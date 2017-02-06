library(dplyr); library(plyr)
library(ggplot2)
library(cowplot)

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
  w = rweibull(np, mu+skew,sd)
  mu-mean(w)+w
}
bimodalDistFunc <- function (n,cpct, mu1, mu2, sig1, sig2) {
  y0 <- rnorm(n,mean=mu1, sd = sig1)
  y1 <- rnorm(n,mean=mu2, sd = sig2)
  
  flag <- rbinom(n,size=1,prob=cpct)
  y <- y0*(1 - flag) + y1*flag 
}
injectBeam <- function(distrib){
  
  np = 1000 # number of particles in bunch
  sddy = 1e-3 #sd of the energy distribution
  w0 = 3; f0 = w0/2/pi # frequency of the reference particle
  p0 = pi/2 #phase of the reference particle
  sdw = w0*sddy
  sdp = p0*2e-2
  
  
  df.p = data.frame(
    wFreq = switch(distrib,
                   "skew" = skewedDistFunc(np, w0, sdw, -2),
                   "bi" = bimodalDistFunc(np,.5,w0-3*sdw,w0+3*sdw,sdw,sdw),
                   "phys" = w0 + 1e3*rnorm(np,sd=sddy)^2 # dw = G*dy^2; w = w0+dw; dy ~ Norm(0,sddy)
    ), 
    Phi = rnorm(np, p0, sdp) # bimodalDistFunc(np,0,p0-3*sdp,p0+3*sdp,sdp,sdp)
  )
  attr(df.p$wFreq, "Synch") <- w0
  attr(df.p$wFreq, "SD") <- sdw
  attr(df.p$Phi, "Synch") <- p0
  attr(df.p$Phi, "SD") <- sdp
  
  df.p
}
Pproj <- function(df, x) colSums(cos(df$wFreq%o%x + df$Phi))

## particle distributions ####
df.p = injectBeam("phys")

.gghist_plot(df.p, "wFreq") -> whist
.gghist_plot(df.p, "Phi") -> phist
plot_grid(whist,phist,nrow=2)

## computing signal ####
w0 = attr(df.p$wFreq, "Synch")
Tstt = 0; Ttot=721; dt = .5/w0 #.5 to satisfy the Nyquist condition
df.s = data.frame(Time=seq(Tstt,Ttot,dt)) %>% mutate(Sgl=Pproj(df.p,Time))

## computing signal peaks ####
p0 = attr(df.p$Phi, "Synch")
Ntot = floor(.5*(w0*Ttot+p0)/pi)
Nstt = floor(.5*(w0*Tstt+p0)/pi)
tnu = (2*pi*Nstt:Ntot-p0)/w0
tnd = tnu+pi/w0

## plotting signal ####
ggplot(df.s, aes(Time, Sgl)) + geom_line(col="gray") + 
  geom_hline(yintercept = c(min(df.s$Sgl), max(df.s$Sgl)), col="red") +
  theme_bw() + 
  geom_point(aes(col=Side), 
             data=data.frame(
               Time=c(tnu,tnd), 
               Sgl=Pproj(df.p,c(tnu,tnd)), 
               Side=rep(c("U","D"),c(length(tnu),length(tnd)))
             ), show.legend = FALSE) +
  scale_color_manual(breaks=c("D","U"), values = c("blue","red")) -> sglplot

## spectral analysis ####
s = ts(df.s$Sgl, start=Tstt, end=Ttot, deltat=dt)
spec.pgram(s,plot = FALSE) -> sps
sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)

x = arrange(sps, desc(Pow))[1:15,]
dw = x[1,"wFreq"]-x[nrow(x),"wFreq"]

sdw = attr(df.p$wFreq,"SD")
ggplot(x,aes(wFreq, Pow))+
  geom_bar(stat="identity", width=dw*1e-2) + 
  theme_bw() + labs(x=expression(omega)) +
  geom_vline(xintercept = w0, col="red") -> fpsplot

## FREQUENCY-SPREAD SIGNAL ERROR HYPOTHESIS ####
# at each point in time the signal is a sum of random variables (b/c wFreq and Phi are RVs)
# and if I recall correctly, the sum of any random variables is a normally distributed
# random variable

test <- function(Ntrl=1000, at){
  p = array(dim=Ntrl)
  for(n in 1:Ntrl) 
    p[n] = Pproj(injectBeam("phys"), at)
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
  print(p)
  # coef
  p
}



test(100, pi/3*100) -> s
gg_qq(s)->tqqplot

library(ggExtra)
ggExtra::ggMarginal(tqqplot, type="density")

## Yup, I was right
if(FALSE){
  tn = (pi)/w0
  ldply(seq(tn-.05,tn+.05,length.out = 10), function(x) {test(100, x)->s; data.frame(At=rep(x,length(s)),Sgl=s)}) -> s
  library(mosaic)
  mean(Sgl~At, data=s) ->xs
  sd(Sgl~At, data=s) -> sds
  cbind(XSgl=xs, SD=sds) %>% as.data.frame() %>% mutate(At = as.numeric(names(xs))) -> s
  
  ggplot(s, aes(x=At, y=XSgl)) + geom_linerange(aes(ymin=XSgl-SD, ymax=XSgl+SD)) + geom_line(col="gray")
}

