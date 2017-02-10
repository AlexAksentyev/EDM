library(dplyr); library(plyr)
library(ggplot2); library(gridExtra)
library(grid)

rm(list=ls(all=TRUE))

source("./RScripts/definitions.R")
source("./RScripts/CSignal.R")

## FUNCTIONS ####
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
compNullX <- function(Ttot, Tstt, df.p, what = "Envelope"){
  w0 = attr(df.p$wFreq, "Synch")
  p0 = attr(df.p$Phi, "Synch")
  
  .dum <- function(Time) floor((w0*Time+p0-pi/2)/2/pi)
  Ntot = .dum(Ttot)
  Nstt = .dum(Tstt)
  d = switch(what, "Envelope" = pi/2, "Node" = 0)
  tnu = (2*pi*Nstt:Ntot-p0+d)/w0; tnu <- tnu[tnu>=0]
  tnd = tnu+pi/w0
  df.pks = data.frame(
    N = c(1:length(tnu),1:length(tnd)),
    Time=c(tnu,tnd), 
    Sgl=Pproj(df.p,c(tnu,tnd)), 
    Side=rep(c("U","D"),c(length(tnu),length(tnd)))
  )
}
compActX <- function(pk.est = df.pks, df.sgl=df.s, what = "Envelope", tol=1e-3){
  require(doParallel); n.cores = detectCores()
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  
  what.f <- eval(parse(text = paste0("which.", switch(what, "Envelope"="max","Node"="min"))))
  
  adply(pk.est, 1, function(s, df.sgl, what.f, tol){
    s[,"Time"] -> x
    filter(df.sgl, Time>x-.5 & Time<x+.5) -> y
    
    dy = y[nrow(y),]-y[1,]
    spline(y$Time, y$Sgl, n=dy$Time/tol) %>% as.data.frame() -> y
    names(y) <- c("Time","Sgl")
    
    slice(y, what.f(abs(Sgl)))
  }, df.sgl, what.f, tol, .parallel = FALSE, .paropts = list(.packages="dplyr")) -> pks
  
  stopCluster(clus)
  
  pks
}

sgl = popPS(CSignal()) ## set up beam phase space

## particle distributions ####
.gghist_plot(sgl@PS, "wFreq") + labs(x=expression(omega)) -> whist
.gghist_plot(sgl@PS, "Phi") + labs(x=expression(phi))-> phist
# grid.arrange(whist,phist)

## computing signal ####
Tstt = 0; Ttot=2000; dt = .25/sgl@wFreq0 # pi/w0 to satisfy the Nyquist condition
sgl <- Signal(sgl, seq(Tstt, Ttot, dt))

## fitting signal ####
f = Sgl ~ nrow(df.p)* exp(lam*Time) * sin(w*Time + p0)
nls(f, data=df.s, start=list(lam=-1.4e-3, w=attr(df.p$wFreq, "Synch"))) -> mod1
mutate(df.s, fSgl = fitted(mod1)) -> df.s

## computing envelopes ####
compNullX(Ttot, Tstt, df.p, "Node")%>%mutate(E="Null")->df.pks0
compActX(df.pks0, what="Node",tol=1e-6) %>%mutate(E="Act")-> df.pks1
df.pks = rbind(df.pks0, df.pks1)

x = df.pks1$Time; x <- c(x[1], x[-length(x)])
df.pks1 %>% mutate(DT = Time - x) %>% filter(N>1)->df.pks1
ggplot(df.pks1) + geom_density(aes(DT))
sd(df.pks1$DT)

## plotting signal ####
ggplot(df.s, aes(Time, Sgl)) + geom_line(col="red",lwd=.05) + 
  theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
  theme(legend.position="top")+
  geom_point(aes(col=E), size=.1, data=df.pks, show.legend = FALSE) +
  scale_color_manual(name="Envelope", values = c("black","blue")) +
  annotation_custom(tableGrob(formatC(coef(summary(mod1))[,1:2],3,format="e")), 
                    xmin=1000, xmax=2500,ymin=-1000, ymax=-650) -> sglplot

t0 = 1900; dt0=25
ggplot(filter(df.s, Time>t0&Time<t0+dt0), aes(Time,Sgl)) + geom_line(col="red") + theme_bw() +
  scale_color_manual(name="Envelope",values=c("black","blue"))+
  geom_point(aes(Time, Sgl, col=E), size=1, data=filter(df.pks, Time>t0&Time<t0+dt0))+labs(x="",y="") -> sglslc

vp <- viewport(width = .4, height = .3, x = .6, y = .5, just = c("right","center"))
# print(sglplot); print(sglslc,vp=vp)

## spectral analysis ####
s = ts(df.s$Sgl, start=Tstt, end=Ttot, deltat=dt)
spec.ar(s,plot = FALSE) -> sps
sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)

x = arrange(sps, desc(Pow))[1:20,]
dw = x[1,"wFreq"]-x[2,"wFreq"]

sdw = attr(df.p$wFreq,"SD")
ggplot(x,aes(wFreq, Pow))+scale_y_continuous(labels=.fancy_scientific) +
  geom_bar(stat="identity", width=dw*.1) + 
  theme_bw() + labs(x=expression(omega)) +
  geom_vline(xintercept = w0, col="red") -> fpsplot

# grid.arrange(sglplot, fpsplot)


## TESTING OUT GROWING PHASE SPACE #### 
if(FALSE){
  ldply(seq(0,721*2,dt), function(x) {
    Pproj(phaseSpace("phys",x),x)->s
    data.frame(Time=x,Sgl=s)
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

## FREQUENCY CREEP ####
if(FALSE){

  names(df.pks0)[c(2:3,5)] <- names(df.pks)[c(2:3,5)]%>%paste0("0")
  names(df.pks1)[c(4:5,3)] <- names(df.pks)[c(2:3,5)]%>%paste0("1")
  join(df.pks0, df.pks1, by=c("N","Side"))%>%dplyr::select(Side,N,E0,Time0,Sgl0,E1,Time1,Sgl1,-fSgl) -> pks
  
  mutate(pks, DT = Time0-Time1, wR = 2/pi*DT/(1+4*N))->pks
  
  lm(wR ~ Time0, data=filter(pks, Time0>500)) -> mod3
  print(coef(summary(mod3)))
  
  pks %>%
    filter(Side=="U") %>% 
    ggplot(aes(Time1, wR)) + geom_line(col="red")+
    theme_bw() + labs(
      y=expression(frac(1,omega[0])-frac(1,omega)),
      x="Time"
    ) #+ 
    # geom_smooth(method="lm", size=.3, col="black") +
    annotation_custom(tableGrob(formatC(coef(summary(mod3))[,1:2], 3, format="e"))) 

}

## INVESTIGATE W ####
if(FALSE){
  library(mosaic)
  mutate(df.s, Part = derivedFactor(
    "<500" = Time < 500,
    "<1000" = Time < 1000,
    "<1500" = Time < 1500,
    "<2000" = Time < 2000,
    .default = "0",
    .method = "first"
  )) -> df.s
  
  library(doParallel); n.cores = detectCores()
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  xpts <- c("f","w0", "df.p","p0")
  clusterExport(clus, xpts)
  
  ddply(df.s, "Part", function(prt) {
    nls(f, data=prt, start = list(w=w0, lam=1e-4))%>%summary%>%coef -> x
    data.frame(w=x["w",1],SEw=x["w",2], lam=x["lam",1], SElam=x["lam",2])
  }, .parallel = FALSE, .paropts = list(.packages="dplyr")) -> wfits
  
  stopCluster(clus)
  
  wfits%>%filter(Part!="0")%>%mutate(dw=w-w0) %>%
    ggplot(aes(Part, dw)) + geom_pointrange(aes(ymin=w-SEw,ymax=w+SEw)) + theme_bw() + labs(y=expression(omega-omega[0])) +
    scale_y_continuous(labels=.fancy_scientific)
  
}
