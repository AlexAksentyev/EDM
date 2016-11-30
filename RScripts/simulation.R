rm(list=ls(all=TRUE))
library(dplyr); library(plyr)
library(reshape2)
fancy <- function(x) formatC(x, 4, format = "e")

## I should show that my formula for var omega is the correct one
## define the measurement model
P = .8; N0 = 6730
errS = 1e-1*N0*P #absolute measurement error

phi = pi/2 #initial signal phase
w = 2*pi*.005#2*pi*120e3*.74e-4 # estimated frequency
.dcs <- function(x) N0 * (1 + P*sin(w*x + phi)) # structural model/expectation

w.guess = rnorm(1, w, .1*w) # if we do modulated sampling, we must guess the signal frequency; this line models the guessing

Ttot = 1000 #total measurement time
fs = w/(2*pi)*sqrt(4.3) #sampling frequency
t = seq(0, Ttot, by = 1/fs);#uniform sampling time 
t2 = seq(0, Ttot, by = .33/fs); plot(.dcs(t2)~t2, type="o", xlab="Time, sec", ylab="Signal, counts"); points(.dcs(t)~t, col="blue")
t2 <- t2[.dcs(t2) > N0*(1 - P*.5) & .dcs(t2) < N0*(1 + P*.5)]; points(.dcs(t2)~t2, col="red"); 
abline(h=N0*(1 + c(-.5,0,.5)*P), col="red", lty=c(2,3,2))
x = seq(0,Ttot,length.out=1000); lines(.dcs(x)~x, lty=4)
min(length(t), length(t2))->l
t <- t[1:l]; t2 = t2[1:l]
t <- list("uni" = t, "mod" = t2)
rm(t2)

Es = llply(t, function(x) {.dcs(x) -> e; names(e) <- x; e}) # expected values for the two samplings


#### modeling ####
Ntrial = 500 #number of trials
matrix(rnorm(Ntrial*l, sd=errS), ncol=Ntrial) -> s
llply(Es, function(e){e+s -> x; rownames(x) <- names(e); x})->s

melt(s) -> s; names(s) <- c("Time","Trl","Sgl", "Type")

## for all trials ####
s %>% ddply(.(Type, Trl), function(trl){
  nls(Sgl~ n*(1 + p*sin(frq*Time + phi)), data=trl, start=list(n = 5000, p = 1, frq=w.guess)) -> m3
  .est = coef(m3); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m3))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}) -> .stats

.se <- function(x) sd(x)/sqrt(length(x))

library(ggplot2)
ggplot(.stats, aes(X.frq, col=Type)) + geom_density() + 
  geom_vline(xintercept=w, linetype=2)+ 
  theme_bw() + labs(x=expression(hat(omega)))

## comparing the efficiencies
daply(.stats, "Type", function(x) .se(x$X.frq)) -> x
x["uni"]/x["mod"]
