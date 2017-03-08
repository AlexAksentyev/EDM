library(dplyr); library(ggplot)
library(reshape2)

rm(list=ls(all=TRUE))

source("./RScripts/RCBunch.R")
source("./RScripts/RCSignal.R")

f0 = .5; w0=f0*2*pi; dt = .37*pi/w0
x = seq(0, 100, by= dt)
e = rnorm(length(x), sd=3e-2)
xs0 = sin(w0*x)
xs1 = sin((w0+1e-9)*x)

s0 = ts(xs0 +e, start=x[1], end=x[length(x)], deltat=dt)
s1 = ts(xs1 + e, start=x[1], end=x[length(x)], deltat=dt)

pspectrum(s0, plot=TRUE) -> fps0
df = data.frame(Pow=fps0$spec, Freq=fps0$freq)
ggplot(df%>%filter(Freq>.45, Freq<.55), aes(Freq, Pow)) + geom_bar(stat="identity")
df[which.max(df$Pow),]
