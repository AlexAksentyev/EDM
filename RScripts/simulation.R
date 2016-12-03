rm(list=ls(all=TRUE))
library(dplyr); library(plyr)
library(reshape2)
library(ggplot2)
library(parallel); library(doParallel); registerDoParallel(detectCores()-1)

fancy <- function(x) formatC(x, 4, format = "e")
.se <- function(x) sd(x)/sqrt(length(x))

.dcs <- function(x, w = w0) N0 * (1 + P*sin(w*x + phi)) # structural model/expectation

.fit <- function(trl){
  nls(Sgl~ n*(1 + p*sin(frq*Time + phi)), data=trl, start=list(n = 5000, p = 1, frq=w.g)) -> m3
  .est = coef(m3); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m3))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}
.simul <- function(data) data %>% ddply(.(Type, Trl), .fit, .parallel=FALSE) %>% 
  mutate(B.frq = X.frq - w0, RB.frq = B.frq/w0, Type = as.factor(Type)) 

.sample <- function(Nprd){
  Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
  assign("Ttot", Ttot, envir = .GlobalEnv)
  t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
  
  ## modulated sampling if the measurements are taken uniformly in time
  t2 = seq(0, Ttot, by = comptn/fs); x = .dcs(t2, w=w.g) #sampled measurements
  t2 <- t2[x > N0*(1 - P*comptn) & x < N0*(1 + P*comptn)] #pick the good points
  rm(x)
  min(length(t1), length(t2))->l
  t1 <- t1[1:l]; t2 = t2[1:l]
  x1 = .dcs(t1); x2 = .dcs(t2)
  names(x1) <- t1; names(x2) <- t2
  list("uniform" = x1, "modulated" = x2)
}
#### preparations ####
## I should show that my formula for var omega is the correct one
w0=3; phi=pi/36; N0=6730; P=.8 # signal
fs=5000; comptn=.33; w.g=rnorm(1,3,0) ## sampling
Es = .sample(3) # expected values for the two samplings
l = length(Es[[1]])
## modeling ##
Ntrial = 50 #number of trials
errS = 3e-2*N0*P #absolute measurement error

## trials differ by error ####
matrix(rnorm(Ntrial*l, sd=errS), ncol=Ntrial) -> s
llply(Es, function(e){e+s -> x; rownames(x) <- names(e); x}) %>% melt ->s; names(s) <- c("Time","Trl","Sgl", "Type")
.stats <- .simul(s)

## plotting results ##
ggplot(.stats, aes(B.frq, col=Type)) + geom_density() +
  theme_bw() + labs(x=expression(bar(omega) - omega))

## comparing the efficiencies
daply(.stats, "Type", function(x) .se(x$X.frq)) -> x
print(x)
c("SE ratio" = x["uniform"]/x["modulated"], "wg/w0" = w.g/w0) %>% print

x = filter(s, Trl==1) %>% mutate(Sgl.o = .dcs(Time), errY = Sgl-Sgl.o) %>% slice(seq(1,2*l, length.out=250))
ggplot(x, aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
  geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, by=1/fs)) %>% mutate(Sgl=.dcs(Time)))

ggplot(.stats) + geom_boxplot(aes(Type, B.frq)) + theme_bw() + labs(y=expression(sigma[hat(omega)]))

## checking the growth of omega se with total time ####
## take the first super test and compute the se for this case
## then exclude the last m periods, recalculate the se
## repeat this reduction-recomputation
## plot results
Es <- .sample(480); length(Es[[1]])->l
rnorm(l, sd=errS)  -> s
ldply(Es, function(e) data.frame("Time"= as.numeric(names(e)), "Sgl"=e+s), .id="Type") %>% mutate(Trl=1)->s
(1:10)*l/10 -> m
llply(m, function(i) s%>%ddply("Type", function(x) slice(x,1:i))) -> s
ldply(s, function(df) .simul(df)%>%mutate(NUM = nrow(df))) %>% 
  mutate(MEAN.SE = SE.frq/sqrt(NUM))->.stats

ggplot(.stats, aes(NUM, MEAN.SE, col=Type)) + geom_point() + 
  scale_y_log10() + scale_x_log10() +
  theme_bw() + labs(x="Sample size",y=expression(sigma[bar(omega)])) + theme(legend.position="top")
.stats %>% group_by(Type)

filter(.stats, Type=="uniform")->x
x[1,c("NUM","MEAN.SE")]/x[9,c("NUM","MEAN.SE")]
write.table(.stats, file="CumFit.txt", quote=FALSE)

## question ##
x = filter(.stats, Type=="uniform") %>% mutate(Time = NUM/fs)
nls(MEAN.SE~I(1/sqrt(Time)), data=x)
