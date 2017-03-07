library(plyr)
library(ggplot2); library(gridExtra)
library(grid)
library(mosaic)

rm(list=ls(all=TRUE))

lblfnt=16
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

source("./RScripts/RCBunch.R")
source("./RScripts/RCSignal.R")

b0 = RCBunch$new(WDist="norm")
b1 = RCBunch$new(WDist="phys")

x0 = diff(quantile(b0$EnsPS[,"wFreq"], c(pnorm(-1), pnorm(1))))
x1 = diff(quantile(b1$EnsPS[,"wFreq"], c(0, pnorm(1)-pnorm(-1))))
x0/x1

dt = .5*pi/b0$Synch["wFreq"]
s0 = RCSignal$new(b0, seq(0, 6000, dt))
s1 = RCSignal$new(b1, seq(0,6000,dt))

w0 = b0$Synch["wFreq"]; p0 = b0$Synch["Phi"]
Ntot = (w0*20000+p0-pi/2)/pi #The number of peaks during a given time
stime <- ((.5+2*(0:Ntot))*pi-p0)/w0 #peak times
stime <- stime[seq(1, length(stime), length.out = 5)] #pick out only several

data.table(WDist="Norm", Time=stime, t(b0$Phase(stime)))%>%melt.data.table(id.vars=c("WDist","Time"))->ps0
data.table(WDist="Chi2", Time=stime, t(b1$Phase(stime)))%>%melt.data.table(id.vars=c("WDist","Time"))->ps1

rbind(ps0, ps1) -> ps; rm(ps0, ps1)

ps[,Ph0 := Time*b0$Synch["wFreq"]+b0$Synch["Phi"]]
ps[,dPh := value-Ph0]
ps[,Shade:=derivedFactor(
  "add" = dPh < -1.5*pi| (dPh > -.5*pi & dPh < .5*pi) | (dPh > 1.5*pi & dPh < 2.5*pi),
  .default="subtract"
)]
ps[, `:=`(S=sin(value),S0=sin(Ph0))]
ps[,`:=`(SS=sum(S)/333, Ph50=median(dPh)), by=c("WDist", "Time")]

ggplot(ps%>%filter(dPh/pi>-2, dPh/pi<3)) + geom_density(aes(dPh/pi)) + 
  geom_density(aes(S), col="red", linetype=2) + 
  facet_grid(Time~WDist, scales = "free_y") +
  labs(x=expression(Delta~Theta/pi)) +
  geom_vline(aes(xintercept=SS), col="darkgreen") +
  thm + 
  geom_segment(aes(x=dPh/pi, xend=dPh/pi, y=0, yend=.1, col=Shade, alpha=.2))


fp = list(func = ValNs ~ 1000 * exp(lam*(Time-d)) * sin(w*Time + p0),
          guess = list(lam=-1.4e-3, w=3, d=100))

s0$fit(fp); #s1$fit(fp)
df = rbind(s0$Signal[,.(Time, ValNs, Fit)][,WDist:="Norm"], s1$Signal[,.(Time, ValNs, Fit)][,WDist:="Chi2"])
df[seq(1, nrow(df), by=2), Fit := df[seq(2,nrow(df), by=2), Fit]] 

ggplot(df) + geom_line(aes(Time, ValNs)) + 
  geom_point(aes(Time, Fit), col="red", size=.5) + 
  facet_grid(WDist~.) + thm + labs(y="Signal")


## CM ####
.gghist <- function(df, name){
  ggplot(df, aes_string(name)) + geom_histogram(aes(y=..density..), col="black",fill="white") + thm +
    geom_density(col="red")
}
b0 = RCBunch$new(SDdy=1e-5)
w0 = b0$Synch["wFreq"]; p0 = b0$Synch["Phi"]
dt = .37*pi/b0$Synch["wFreq"]
s0 = RCSignal$new(b0, seq(0, 3000, dt))

## phase space 
df = data.table(wFreq=b0$EnsPS[,"wFreq"], Phi=b0$EnsPS[,"Phi"])
whist <- .gghist(df, "wFreq") + labs(x=expression(omega))
phist <- .gghist(df, "Phi") + labs(x=expression(phi))
grid.arrange(whist, phist, nrow=2)

## signal and fit

fp = list(func = ValNs ~ 1000 * exp(lam*(Time-d)) * sin(w*Time + p0),
          guess = list(lam=-1.4e-3, w=3, d=100))

s0$fit(fp) -> mod3l
s0$findPts("Envelope", w.guess=s0$ModelCoeff["w",1])
df = s0$specPts[Which=="Optim",.(N, Time, Side, Val)]
df[,`:=`(Fit=predict(mod3l, newdata=list("Time"=Time)))]

ggplot(s0$Signal[seq(1,nrow(s0$Signal), length.out=3102)]) + 
  geom_line(aes(Time, Val), size=.15, col="red") +
  geom_point(aes(Time, Val), data=df, size=.25) +
  geom_point(aes(Time, Fit), data=df, size=.15, col="blue") +
  thm + labs(y="Signal")

s0$Spectrum()
