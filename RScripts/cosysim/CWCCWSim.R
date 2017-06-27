library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)

Nrev = 25000 #simulation revolutions
Nrps = 1e6 #revolutions per second
sigtilt = 1e-3 #increment

.complete <- function(DT){
  DT[,':='(
    ThXY = atan(Sy/Sx),
    ThXZ = atan(Sz/Sx),
    ThYZ = atan(Sy/Sz)
  )
  ]
  DT[,':='(
    Wz = ThXY/Nrev*Nrps,
    Wy = ThXZ/Nrev*Nrps,
    Wx = ThYZ/Nrev*Nrps
  )]
  
}

.prep <- function(Direct, Reverse){
  .complete(Direct)
  .complete(Reverse)
  
  names(Direct) <- paste0(names(Direct),".CW")
  names(Reverse) <- paste0(names(Reverse),".CCW")
  
  Direct[rep(1:nrow(Direct),each=nrow(Reverse)),]->Direct
  Reverse[rep(1:nrow(Reverse),nrow(Reverse)),]->Reverse
  cbind(Direct, Reverse) -> DR
  DR[,':='(
    dWx = Wx.CW+Wx.CCW,
    dWy = Wy.CW-Wy.CCW,
    dWz = Wz.CW+Wz.CCW
  )]
  
  return(DR)
}

ldply(1:5, function(k){
  Direct <- data.table(read_csv(paste0("~/git/COSYINF/test/fort.",(6+2*k))))
  Reverse <- data.table(read_csv(paste0("~/git/COSYINF/test/fort.", (6+2*k+1))))
  
  DR <- .prep(Direct, Reverse)
    
  DR[,sigTILT:=(k-1)*sigtilt]
  DR
}) %>% data.table() -> DR

DR[,sigTILT:=as.factor(DR$sigTILT)]

##X-DGAMMA OMEGA PLOT 
dr <- DR[,.(sigTILT, dgamma.CW, dgamma.CCW, x.CW, x.CCW, Wy.CW, Wy.CCW, Wx.CW, Wx.CCW)]
wch = "CCW"
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wy.",wch), col=paste0("dgamma.",wch), shape="sigTILT")) + 
  theme_bw() + labs(x="x",y=expression(Omega[y]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p1
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wx.",wch), col=paste0("dgamma.",wch), shape="sigTILT")) + 
  theme_bw() + labs(x="x",y=expression(Omega[x]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p2
plot_grid(p1,p2,nrow=2)

## HISTOGRAMS
ggplot(dr) + geom_histogram(aes(Wx.CW, fill=dgamma.CW),binwidth=1e-4) + facet_grid(sigTILT~.)

##CROSS PLOT
filter(DR, abs(dWy)>0)%>% 
  ggplot(aes(abs(dWy),abs(dWx), col=sigTILT))+
  geom_point() + theme_bw() + theme(legend.position="top") +
  labs(x=expression(abs(Delta~Omega[y])), y=expression(abs(Sigma~Omega[x])), col=expression(sigma[tilt])) +
  geom_smooth(method="lm", se=FALSE)

ddply(DR, .(sigTILT), function(df){
  lm(abs(dWx)~abs(dWy), df, abs(dWy)>0) -> m
  coef(summary(m))[1:2, 1:2]-> mt
  c(Itct=mt[1,1],SEItct=mt[1,2],Slp=mt[2,1],SESlp=mt[2,2])
}) -> .stats

ggplot(.stats, aes(sigTILT, Itct)) + geom_point() +
  geom_errorbar(aes(ymin=Itct-SEItct, ymax=Itct+SEItct)) +
  theme_bw() + labs(x=expression(sigma[tilt]), y=expression(abs(Sigma~Omega[x])))
