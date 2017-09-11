library(readr)
library(lattice)
library(ggplot2)
library(data.table)
library(dplyr)
library(mosaic)

## CONSTANTS ####
clight = 299792458 #m/s
Lcir = 145.85 #m
e=1.602e-19 #C
md = 2*1.672e-27 #kg
G=-.14

## INPUT ####

filename = "~/git/COSYINF/test/fort.10"
TSSdata <- data.table(read_csv(filename, skip=2))
para <- read_csv(filename, n_max=1)
setattr(TSSdata,"Parameter", para)

Ntps = clight*para$beta/Lcir
alpha <- e/md/para$gamma * (para$gamma*G+1)


TSSdata[,`:=`(W=2*pi*tune*Ntps, Dim=para$Dimension)]

TSSdata[,Wx:=W*nx]

## PLOTS ####

D0 <- TSSdata[y==1e-8&NRG==1e-8][,.(meanWx=mean(Wx)),by=c("muTILT","x")]
D1 <- TSSdata[x==1e-8&y==1e-8&NRG==1e-8]
D2 <- TSSdata[y==1e-8&muTILT==0][,.(meanWx=mean(Wx)),by=c("NRG","x")]

wireframe(meanWx~muTILT*x, data=D0,drape=TRUE)
wireframe(meanWx~NRG*x, data=D2,drape=TRUE)

ddply(D0,"x", function(df){
  coef(summary(lm(meanWx~muTILT, data=df)))[2,1:2]
}) -> m4c

ggplot(m4c, aes(x, Estimate)) + geom_point() +
  # geom_pointrange(aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`)) +
  labs(y=expression(Delta~bar(Omega)[x]/Delta~bar(Theta))) +
  theme_bw()


ggplot(D0[x==1e-8], aes(muTILT,meanWx)) + geom_point() + theme_bw()
ggplot(D0, aes(x,meanWx)) + geom_point(aes(col=muTILT))

ggplot(D1[muTILT==0]) + geom_histogram(aes(Wx),fill="white",col="black") 

