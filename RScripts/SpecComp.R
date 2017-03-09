library(dplyr); library(ggplot2)
library(reshape2)
library(psd)

rm(list=ls(all=TRUE))

source("./RScripts/RCBunch.R")
source("./RScripts/RCSignal.R")

f0 = .5; w0=f0*2*pi; dt = .37*pi/w0
x = seq(0, 1080, by= dt)

registerDoParallel(detectCores()-1)
ldply(c(c(1,2,5)%o%10^(-9:-2)), function(dw, x){
  e = rnorm(length(x), sd=3e-2)
  xs0 = sin((w0-dw)*x)
  xs1 = sin((w0+dw)*x)
  
  s0 = ts(xs0 +e, start=x[1], end=x[length(x)], deltat=dt)
  s1 = ts(xs1 + e, start=x[1], end=x[length(x)], deltat=dt)
  
  bspec(s0)->ps0
  bspec(s1)->ps1
  
  q0 <- quantile(ps0, probs = seq(.1,.9,.1))
  q1 <- quantile(ps1, probs = seq(.1,.9,.1))
  
  ldply(1:9, function(x) (lm(q1[,x]~q0[,x])%>%summary%>%coef)[1,1:2]) -> stats
  
  transmute(stats, DW=dw, Q=seq(.1,.9,.1), Estimate=Estimate, SE=`Std. Error`)
}, x, .parallel=TRUE) -> stats

ggplot(stats%>%filter(), aes(DW, Estimate)) + 
  geom_pointrange(aes(ymin=Estimate-SE, ymax=Estimate+SE)) +
  geom_point(aes(DW,DW), col="red") + 
  scale_x_log10() + scale_y_log10() + 
  theme_bw() + labs(x=expression(Delta~omega)) +
  theme(legend.position="top")
