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

df = data.frame(wFreq=w, Phi=p)
df <- data.frame(Time=seq(0,500,.5/w0)) %>% mutate(Sgl=Pproj(df,Time))

ggplot(df, aes(Time, Sgl)) + geom_line() + geom_hline(yintercept = c(min(df$Sgl), max(df$Sgl)), col="red") +
  theme_bw()

