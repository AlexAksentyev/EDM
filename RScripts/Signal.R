library(dplyr); library(plyr)
library(ggplot2)


np = 1000 # number of particles in bunch
w0 = 3 # frequency of the reference particle
p0 = pi/3 #phase of the reference particle
w = rweibull(np, 2, w0*3e-3) # particle spin precession frequencies
w <- w0-mean(w) + w # to center the distribution on w0 (actually centers only when shape=scale = Gauss)
hist(w, freq=FALSE); lines(density(w,adjust = 2)); lines(density(rnorm(np, w0, sd(w)),adjust=2),col="red", lty=2)
p = rnorm(np, p0, 1e0) # initial phases

Pproj <- function(df, x) colSums(cos(df$wFreq%o%x + df$Phi))

Tstt = 0; Ttot=500
df.p = data.frame(wFreq=w, Phi=p)
df.s = data.frame(Time=seq(Tstt,Ttot,.5/w0)) %>% mutate(Sgl=Pproj(df.p,Time))

Ntot = floor(.5*(w0*Ttot+p0)/pi)
Nstt = floor(.5*(w0*Tstt+p0)/pi)
tn = (2*pi*Nstt:Ntot-p0)/w0

ggplot(df.s, aes(Time, Sgl)) + geom_line() + geom_hline(yintercept = c(min(df.s$Sgl), max(df.s$Sgl)), col="red") +
  theme_bw() + geom_point(aes(col="red"), data=data.frame(Time=tn, Sgl=Pproj(df.p,tn)), show.legend = FALSE)


# the signal frequency is determined by the MODE of the distribution