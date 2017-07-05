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

TSSdata[,`:=`(Omega=mu*Ntps, Dim=para$Dimension)]


alpha = e/md/clight*(para$gamma*G+1)/para$gamma


## PLOTS ####

TSSdata[x==0&y==0&ERG==min(abs(ERG))]%>%ggplot(aes(muTILT, nx))+
  geom_point()+
  theme_bw()
