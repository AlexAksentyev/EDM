library(dplyr); library(plyr)
library(ggplot2); library(gridExtra)
library(grid)

source("./RScripts/definitions.R")

.gghist_plot <- function(df, name){
  ggplot(df, aes_string(name)) +
    geom_histogram(aes(y=..density..), fill="white", color="black") +
    geom_density() + 
    stat_function(fun=dnorm, 
                  args=list(mean=attr(df[,name],"Synch"), 
                            sd=attr(df[,name],"SD")), 
                  col="red", lty=2) +
    theme_bw()
}
skewedDistFunc <- function(n, mu, sd, skew){
  w = rweibull(n, mu+skew,sd)
  mu-mean(w)+w
}
bimodalDistFunc <- function (n,cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n,mean=mu1, sd = sig1)
  y1 <- rnorm(n,mean=mu2, sd = sig2)
  
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag 
}
phaseSpace <- function(distrib, at=0, np=1000){
  
  sddy = 1e-3 + 5e-6*at #sd of the energy distribution
  G = 1e3 # dw = G*dy^2
  w0 = 3 # frequency of the reference particle
  p0 = 0 #phase of the reference particle
  sdw = G*sddy^2 #for all distributions other than phys
  sdp = 2e-2
  
  
  df.p = data.frame(
    wFreq = switch(distrib,
                   "norm" = rnorm(np, w0, sdw),
                   "skew" = skewedDistFunc(np, w0, sdw, 2),
                   "bi" = bimodalDistFunc(np,.5,w0-3*sdw,w0+3*sdw,sdw,sdw),
                   "phys" = w0 + G*rnorm(np,sd=sddy)^2 # dw = G*dy^2; w = w0+dw; dy ~ Norm(0,sddy)
    ), 
    Phi = rnorm(np, p0, sdp) # bimodalDistFunc(np,0,p0-3*sdp,p0+3*sdp,sdp,sdp)
  )
  attr(df.p$wFreq, "Synch") <- w0
  attr(df.p$wFreq, "SD") <- sdw
  attr(df.p$Phi, "Synch") <- p0
  attr(df.p$Phi, "SD") <- sdp
  
  df.p
}
Pproj <- function(df, x) colSums(sin(df$wFreq%o%x + df$Phi))

## particle distributions ####
df.p = phaseSpace("norm")

.gghist_plot(df.p, "wFreq") + labs(x=expression(omega)) -> whist
.gghist_plot(df.p, "Phi") + labs(x=expression(phi))-> phist
grid.arrange(whist,phist)

## computing signal ####
w0 = attr(df.p$wFreq, "Synch")
Tstt = 0; Ttot=721*3; dt = .5/w0 #.5 to satisfy the Nyquist condition
df.s = data.frame(Time=seq(Tstt,Ttot,dt)) %>% mutate(Sgl=Pproj(df.p,Time))

## fitting signal ####
p0 = attr(df.p$Phi, "Synch")
f = Sgl ~ nrow(df.p)* exp(lam*Time) * sin(w*Time + p0)
nls(f, data=df.s, start=list(lam=-1.4e-3, w=w0)) -> mod1; print(coef(summary(mod1)))
mutate(df.s, fSgl = fitted(mod1)) -> df.s

## computing signal peaks ####
dum <- function(Time) floor((w0*Time+p0-pi/2)/2/pi)
Ntot = dum(Ttot)
Nstt = dum(Tstt)
tnu = (2*pi*Nstt:Ntot-p0+pi/2)/w0
tnd = tnu+pi/w0
df.pks = data.frame(
  Time=c(tnu,tnd), 
  Sgl=Pproj(df.p,c(tnu,tnd)), 
  Side=rep(c("U","D"),c(length(tnu),length(tnd)))
)

## plotting signal ####
ggplot(df.s, aes(Time, Sgl)) + geom_line(col="gray",lwd=.05) + 
  theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
  geom_point(size=.1,data=df.pks, show.legend = FALSE) +
  # scale_color_manual(breaks=c("D","U"), values = c("blue","red")) +
  # geom_point(aes(Time, fSgl), size=.5, shape="x", col="magenta") +
  annotation_custom(tableGrob(formatC(coef(summary(mod1))[,1:2],3,format="e")), 
                    xmin=1000, xmax=2500,ymin=-1000, ymax=-650) -> sglplot

ggplot(filter(df.s, Time>1000&Time<1010), aes(Time,Sgl)) + geom_line(col="gray") + theme_bw() +
  geom_point(aes(Time, Sgl), size=1, data=filter(df.pks, Time>1000&Time<1010))+labs(x="",y="") -> sglslc

vp <- viewport(width = .4, height = .3, x = .6, y = .5, just = c("right","center"))
print(sglplot); print(sglslc,vp=vp)

## spectral analysis ####
s = ts(df.s$Sgl, start=Tstt, end=Ttot, deltat=dt)
spec.pgram(s,plot = FALSE) -> sps
sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)

x = arrange(sps, desc(Pow))[1:15,]
dw = x[1,"wFreq"]-x[2,"wFreq"]

sdw = attr(df.p$wFreq,"SD")
ggplot(x,aes(wFreq, Pow))+scale_y_continuous(labels=.fancy_scientific) +
  geom_bar(stat="identity", width=dw*.1) + 
  theme_bw() + labs(x=expression(omega)) +
  geom_vline(xintercept = w0, col="red") -> fpsplot

plot_grid(sglplot,fpsplot,nrow=2)

## FREQUENCY-SPREAD SIGNAL ERROR HYPOTHESIS ####
# at each point in time the signal is a sum of random variables (b/c wFreq and Phi are RVs)
# and if I recall correctly, the sum of any random variables is a normally distributed
# random variable

test <- function(Ntrl=1000, at){
  p = array(dim=Ntrl)
  for(n in 1:Ntrl) 
    p[n] = Pproj(phaseSpace("phys", at), at)
  p
}
gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
                  labels = names(x)){
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  if(!is.null(labels)){ 
    df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
  
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2], col="red") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
    theme_bw()
  
  if(!is.null(labels)) p <- p + geom_text( aes(label = label))
  # print(p)
  # coef
  p
}



test(100, (pi+pi/3)/w0*100) -> s
gg_qq(s)->tqqplot

library(ggExtra)
ggExtra::ggMarginal(tqqplot, type="density")

## Yup, I was right. Or is it because the distribution of phi is normal? 
# just looking at the error bars from phaseSpace
if(FALSE){ 
  tn=pi/w0
  ldply(seq(tn-5e-3,tn+5e-3,length.out=10), function(x) {test(100,x) ->s; data.frame(At=rep(x,length(s)),Sgl=s)}) -> s
  
  library(mosaic)
  mean(Sgl~At, data=s) -> xs
  sd(Sgl~At, data=s) -> sds
  data.frame(XSgl=xs,SD=sds,At=as.numeric(names(xs))) -> s
  
  ggplot(s, aes(x=At, y=XSgl)) + geom_linerange(aes(ymin=XSgl-SD,ymax=XSgl+SD)) + geom_line(col="gray")+ theme_bw()
}

## TESTING OUT GROWING PHASE SPACE #### 
if(FALSE){
  ldply(seq(0,721*2,dt), function(x) {
    Pproj(phaseSpace("phys",x),x)->s
    data.frame(Time=rep(x,length(s)),Sgl=s)
  }) -> s

  ## i'll try fitting now
  f = Sgl ~ nrow(df.p)* exp(lam*(Time))* sin(w*Time + g*Time^2 + p0)
  nls(f, data=s, start = list(lam=log(.25)/500, w=rnorm(1,w0,1e-4), g=5e-6)) -> mod
  print(summary(mod))
  
  mutate(s, fSgl = fitted(mod)) -> s
  
  ggplot(s, aes(x=Time, y=Sgl)) + geom_line(col="gray",lwd=.05)+ theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
    # geom_point(aes(Time,fSgl), col="magenta", size=.1) +
    annotation_custom(tableGrob(formatC(coef(summary(mod))[,1:2],3,format="e")),
                      xmin = 600, xmax=1500, ymin=-1000, ymax=-500)
  
}



