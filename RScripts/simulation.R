rm(list=ls(all=TRUE))
library(dplyr); library(plyr)
library(reshape2)
library(ggplot2)
library(parallel); library(doParallel); registerDoParallel(detectCores()-1)

## definitions ####
fancy <- function(x) formatC(x, 4, format = "e")
.se <- function(x) sd(x)/sqrt(length(x))

.dcs <- function(x, w = w0) N0 * (1 + P*exp(lam.decoh*x)*sin(w*x + phi)) # structural model/expectation
.ddcs <- function(x, w = w0) N0*P*exp(lam.decoh*x)*cos(w*x + phi) # signal derivative

.compVarF <- function(df, err = errS){
  ftr <- sum(df$Drvt^2)
  mutate(df, Wt = Drvt^2/ftr, WtT = Time*Wt, WtTT = Time^2*Wt)->df
  (sum(df$WtTT) - sum(df$WtT)^2)*ftr -> denom
  err/sqrt(denom)
}

.sample <- function(Nprd, Ntrl = 1){
  Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
  assign("Ttot", Ttot, envir = .GlobalEnv)
  t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
  
  ## modulated sampling if the measurements are taken uniformly in time
  t2 = seq(0, Ttot, by = comptn/fs); x = .dcs(t2, w=w.g) #sampled measurements
  t2 <- t2[x > N0*(1 - P*comptn) & x < N0*(1 + P*comptn)] #pick the good points
  rm(x)
  min(length(t1), length(t2))->l
  t1 <- t1[1:l]; t2 = t2[1:l]
  df1 = data.frame("Type" = "uniform","Time" = t1, "XSgl" = .dcs(t1), "Drvt" = .ddcs(t1))
  df2 = data.frame("Type" = "modulated","Time" = t2, "XSgl" = .dcs(t2), "Drvt" = .ddcs(t2))
  rbind(df1,df2)[rep(seq(2*l),Ntrl),] %>% arrange(Type) %>%
    mutate(Trl = rep(seq(Ntrl), each=l)%>%rep(2), Sgl = XSgl + rnorm(Ntrl*l, sd=errS)%>%rep(2))
}
.fit <- function(trl){
  nls(Sgl~ n*(1 + p*sin(frq*Time + phi)), data=trl, start=list(n = 5000, p = 1, frq=w.g)) -> m3
  .est = coef(m3); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m3))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}

#### preparations ####
## I should show that my formula for var omega is the correct one
w0=3; phi=pi/36; N0=6730; P=.4 ; lam.decoh = log(.1/.4)/1000# signal
fs=5000; comptn=.33; w.g=rnorm(1,w0,0*w0) ## sampling
errS = 3e-2*N0*P #absolute measurement error

## modeling ##

## 1) trials with different error ####
s = .sample(5,15) # sample signal for the two samplings
# for each type of sampling, and each trial, fit
s %>% ddply(.(Type, Trl), function(trl) .fit(trl)%>%mutate(SEAN.frq = .compVarF(trl))) -> .stats

.stats%>%group_by(Type) %>% dplyr::summarise(NUM = n(), SE = sd(X.frq), SEAN = SEAN.frq[1], RD = (SE-SEAN)/SEAN)

## comparing the efficiencies
daply(.stats, "Type", function(x) .se(x$X.frq)) -> x
print(x)
c("SE ratio" = x["uniform"]/x["modulated"], "wg/w0" = w.g/w0) %>% print

## plot of sample points on signal
x = filter(s, Trl==1) 
x %>% mutate(errY = Sgl-XSgl) %>% slice(seq(1,nrow(x), length.out=250)) %>% 
  ggplot(aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
  geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, by=1/fs)) %>% mutate(Sgl=.dcs(Time)))

## 2) checking the growth of omega se with total time ####
s <- .sample(15); nrow(s)/2->l #divide by 2 b/c 2 sampling types

(1:10)*l/10 -> m # split the super-sample into a sequence of cumulatively increasing subsamples;
                  # m contains the lengths of the subsamples 
# for a given length i, pick the first i measurements from each sampling; this is the subsample; combine the samplings
llply(m, function(i) s%>%ddply("Type", function(x) slice(x,1:i))) -> s
# fit the subsamples
ldply(s, function(df) df%>%ddply("Type", .fit)%>%mutate(NUM = nrow(df), SEAN.frq = .compVarF(df))) %>% 
  mutate(SE.MEAN = SE.frq/sqrt(NUM))->.stats

ggplot(.stats, aes(NUM, SE.MEAN, col=Type)) + geom_point() + 
  scale_y_log10() +
  theme_bw() + labs(x="Sample size",y=expression(sigma[bar(omega)])) + theme(legend.position="top")

ggplot(mutate(.stats, Dev = (SE.frq-SEAN.frq)/SEAN.frq), aes(NUM, abs(Dev), col=Type)) + geom_point() + 
  theme_bw() + labs(x="Sample size", y="(Simulation - Analytic)/Analytic SE") + theme(legend.position="top")

ggplot(.stats, aes(SEAN.frq, SE.frq, col=Type)) + geom_point() + 
  theme_bw() + labs(x="Analytic SE", y="Simulation SE") + theme(legend.position="top") + 
  geom_abline(slope=1,intercept=0) + 
  scale_y_log10() + scale_x_log10()

## 3) checking the analytical formula ####
## the formula is
## SEw = SE error / sqrt(sum xj * (sum ti^2 xi/sum xj) - (sum ti xi/sum xj)^2)
## this doesn't account for damping
## !!!!!!! have to see if damping matters, though !!!!!!
s = .sample(480); l = nrow(s)
.stats <- ddply(s,"Type", .fit) %>% mutate(SEAN.frq = daply(s, "Type", .compVarF))

.stats[1,"SE.frq"]/.stats[2,"SE.frq"]
.stats[1,"SEAN.frq"]/.stats[2,"SEAN.frq"]

x = mutate(s, errY = Sgl-XSgl) %>% slice(seq(1,l, length.out=250))
ggplot(x, aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
  geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, length.out=250)) %>% mutate(Sgl=.dcs(Time))) +
  theme(legend.position="top")

 