rm(list=ls(all=TRUE))
library(dplyr); library(plyr)
library(reshape2)
library(ggplot2)
library(parallel); library(doParallel); registerDoParallel(detectCores()-1)

fancy <- function(x) formatC(x, 4, format = "e")
.se <- function(x) sd(x)/sqrt(length(x))

.fit <- function(trl){
  nls(Sgl~ n*(1 + p*sin(frq*Time + phi)), data=trl, start=list(n = 5000, p = 1, frq=w.g)) -> m3
  .est = coef(m3); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m3))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}
.simul <- function(data) data %>% ddply(.(Type, Trl), .fit, .parallel=TRUE) %>% 
  mutate(B.frq = X.frq - w0, RB.frq = B.frq/w0, Type = as.factor(Type)) 

#### preparations ####
## I should show that my formula for var omega is the correct one
## define the measurement model
P = .8; N0 = 6730
errS = 3e-2*N0*P #absolute measurement error

phi = pi/36 #initial signal phase
w0 = 3 #signal frequency
.dcs <- function(x, w = w0) N0 * (1 + P*sin(w*x + phi)) # structural model/expectation

w.g = rnorm(1, w0, 0*w0) # if we do modulated sampling, we must guess the signal frequency; this line models the guessing

Ttot = (3*2*pi-phi)/w0 #total measurement time
fs = 5000 #sampling frequency
t = seq(0, Ttot, by = 1/fs);#uniform sampling time 

## modulated sampling if the measurements are taken uniformly in time
t2 = seq(0, Ttot, by = .05/fs); 
t2 <- t2[.dcs(t2, w = w.g) > N0*(1 - P*.05) & .dcs(t2, w = w.g) < N0*(1 + P*.05)];
min(length(t), length(t2))->l
t <- t[1:l]; t2 = t2[1:l]
t <- list("uniform" = t, "modulated" = t2)
rm(t2)

## If the modulated measurements aren't taken uniformly, but instead are distributed normally
## about the points where the sine is 0, then w*t_n + phi = pi*n => the modulated samling measurement time model is 
## t_n = pi/w * n + phi/w

Es = llply(t, function(x) {.dcs(x) -> e; names(e) <- x; e}) # expected values for the two samplings

## modeling ##
Ntrial = 500 #number of trials

## trials differ by error ####
matrix(rnorm(Ntrial*l, sd=errS), ncol=Ntrial) -> s
llply(Es, function(e){e+s -> x; rownames(x) <- names(e); x}) %>% melt ->s; names(s) <- c("Time","Trl","Sgl", "Type")
.stats <- .simul(s)

## plotting results ##
ggplot(.stats, aes(X.frq, col=Type)) + geom_density() + 
  geom_vline(xintercept=w0, linetype=2) +
  theme_bw() + labs(x=expression(hat(omega)))

## comparing the efficiencies
daply(.stats, "Type", function(x) .se(x$X.frq)) -> x
print(x)
c("SE ratio" = x["uniform"]/x["modulated"], "wg/w0" = w.g/w0) %>% print

ggplot(.stats, aes(B.frq, col=Type)) + geom_density() +
  theme_bw() + labs(x=expression(bar(omega) - omega))

x = filter(s, Trl==1) %>% mutate(Sgl.o = .dcs(Time), errY = Sgl-Sgl.o) %>% slice(seq(1,2*l, length.out=250))
ggplot(x, aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
  geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, by=1/fs)) %>% mutate(Sgl=.dcs(Time)))

ggplot(.stats) + geom_boxplot(aes(Type, B.frq)) + theme_bw() + labs(y=expression(sigma[hat(omega)]))
