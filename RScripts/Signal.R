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
skewedDistFunc <- function(n, mu, sd, skew){
  w = rweibull(np, mu+skew,sd)
  mu-mean(w)+w
}
bimodalDistFunc <- function (n,cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n,mean=mu1, sd = sig1)
  y1 <- rnorm(n,mean=mu2, sd = sig2)
  
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag 
}

np = 1000 # number of particles in bunch
w0 = 3; f0 = w0/2/pi # frequency of the reference particle
p0 = pi/3 #phase of the reference particle
sdw = w0*1e-3
sdp = p0*3e-2
sddy = 1e-3 #sd of the energy distribution
dis = "phys"

## particle distributions ####
df.p = data.frame(
  wFreq = switch(dis,
    "skew" = skewedDistFunc(np, w0, sdw, -2),
    "bi" = bimodalDistFunc(np,.5,w0-3*sdw,w0+3*sdw,sdw,sdw),
    "phys" = w0 + 1e3*rnorm(np,sd=sddy)^2 # dw = G*dy^2; w = w0+dw; dy ~ Norm(0,sddy)
  ), 
  Phi = rnorm(np, p0, sdp) # bimodalDistFunc(np,0,p0-3*sdp,p0+3*sdp,sdp,sdp)
)

par(mfrow=c(2,1))
.hist_plot(df.p$wFreq,w0,expression(omega))
.hist_plot(df.p$Phi,p0,expression(phi))
par(mfrow=c(1,1))

## computing signal ####
Pproj <- function(df, x) colSums(cos(df$wFreq%o%x + df$Phi))

Tstt = 0; Ttot=1000; dt = .5/w0 #.5 to satisfy the Nyquist condition
df.s = data.frame(Time=seq(Tstt,Ttot,dt)) %>% mutate(Sgl=Pproj(df.p,Time))

## computing signal peaks ####
Ntot = floor(.5*(w0*Ttot+p0)/pi)
Nstt = floor(.5*(w0*Tstt+p0)/pi)
tnu = (2*pi*Nstt:Ntot-p0)/w0
tnd = tnu+pi/w0

## plotting signal ####
ggplot(df.s, aes(Time, Sgl)) + geom_line() + 
  geom_hline(yintercept = c(min(df.s$Sgl), max(df.s$Sgl)), col="red") +
  theme_bw() + 
  geom_point(aes(col=Side), 
             data=data.frame(
               Time=c(tnu,tnd), 
               Sgl=Pproj(df.p,c(tnu,tnd)), 
               Side=rep(c("U","D"),c(length(tnu),length(tnd)))
             ), show.legend = FALSE) +
  scale_color_manual(breaks=c("D","U"), values = c("blue","red"))

## spectral analysis ####
s = ts(df.s$Sgl, start=Tstt, end=Ttot, deltat=dt)
spec.pgram(s,plot = FALSE) -> sps
sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)

filter(sps, wFreq>w0-5*sdw, wFreq<w0+5*sdw) %>% ggplot(aes(wFreq, Pow))+
  geom_bar(stat="identity", width=1e-3) + 
  theme_bw() + labs(x=expression(omega)) +
  geom_vline(xintercept = w0, col="red")
