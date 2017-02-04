library(dplyr); library(plyr)
library(ggplot2)

.hist_plot <- function(x, x0, name){
  require(quantmod)
  hist(x, freq=FALSE, main = paste("Histogram of ", name), xlab=name)
  dx = density(x,adjust = 2)
  xm = dx$x[findPeaks(dx$y)[1]]
  lines(dx)
  lines(density(rnorm(length(x), x0, sd(x)),adjust=2),col="red", lty=2)
  legend("topright",lty=c(1,2),col=c("black","red"), legend = c("density","Gauss"))
  # abline(v=c(xm,x0), lty=c(1,2), col=c("black","red"))
  # legend("top", bty="n",
  #        legend=bquote(frac(x[b]-x[s],x[s])~"="~.(formatC((xm-x0)/x0*100,2,format="e"))~"%")
  # )
}

np = 1000 # number of particles in bunch
w0 = 3; f0 = w0/2/pi # frequency of the reference particle
p0 = pi/3 #phase of the reference particle
sdw = w0*3e-3
sdp = p0*3e-1

## particle distributions ####
df.p = data.frame(wFreq = rweibull(np, w0+5,sdw)) %>% # particle spin precession frequencies
  mutate(wFreq = w0-mean(wFreq)+wFreq) %>% # centering on w0 (actually centers only when shape=scale = Gauss)
  mutate(Phi = rnorm(np, p0, sdp)) # initial phases

par(mfrow=c(2,1))
.hist_plot(df.p$wFreq,w0,expression(omega))
.hist_plot(df.p$Phi,p0,expression(phi))
par(mfrow=c(1,1))

## computing signal ####
Pproj <- function(df, x) colSums(cos(df$wFreq%o%x + df$Phi))

Tstt = 0; Ttot=500; dt = .5/w0
df.s = data.frame(Time=seq(Tstt,Ttot,dt)) %>% mutate(Sgl=Pproj(df.p,Time))

## computing peaks ####
Ntot = floor(.5*(w0*Ttot+p0)/pi)
Nstt = floor(.5*(w0*Tstt+p0)/pi)
tn = (pi*Nstt:(2*Ntot)-p0)/w0

## plotting signal ####
ggplot(df.s, aes(Time, Sgl)) + geom_line() + geom_hline(yintercept = c(min(df.s$Sgl), max(df.s$Sgl)), col="red") +
  theme_bw() + geom_point(aes(col="red"), data=data.frame(Time=tn, Sgl=Pproj(df.p,tn)), show.legend = FALSE)

## spectral analysis ####
s = ts(df.s$Sgl, start=Tstt, end=Ttot, deltat=dt)
spec.ar(s) -> sps
sps$freq[which.max(sps$spec)] -> fsgl
abline(v=c(fsgl,f0), col=c("black","red"),lty=2)
legend("topright",lty=2,col=c("black","red"),legend = c("Bunch","Synchronous"))
legend("top", bty="n", 
       legend=bquote(frac(f[b]-f[s],f[s])~"="~.(round((fsgl-f0)/f0*100,2))~"%")
)
