library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)


.complete <- function(DT){
  Nrev = attr(DT,"Parameters")$Nrev #revolutions per simulation
  Nrps = 1e6 #revolutions per second
  
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

ldply(1:6, function(k){
  ddf = paste0("~/git/COSYINF/test/fort.",(6+2*k))
  rdf = paste0("~/git/COSYINF/test/fort.",(6+2*k+1))
  
  para = read_csv(ddf, n_max=1)
  Direct <- data.table(read_csv(ddf, skip=2))
  attr(Direct, "Parameters")<-para
  para =  read_csv(rdf, n_max=1)
  Reverse <- data.table(read_csv(rdf, skip=2))
  attr(Reverse, "Parameters")<-para
  
  DR <- .prep(Direct, Reverse)
    
  DR[,sigTILT:=para$sigTILT]
  DR
}) %>% data.table() -> DR

fsigTILTlevels = unique(DR$sigTILT)
fsigTILTlevels[order(fsigTILTlevels,decreasing=FALSE)]%>%formatC(2,format="e") -> fsigTILTlevels

DR[,fsigTILT:=factor(DR$sigTILT,labels=fsigTILTlevels)]

##X-DGAMMA OMEGA PLOT 
dr <- DR[,.(fsigTILT,sigTILT, dgamma.CW, dgamma.CCW, x.CW, x.CCW, Wy.CW, Wy.CCW, Wx.CW, Wx.CCW)]
wch = "CCW"
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wy.",wch), col=paste0("dgamma.",wch), shape="fsigTILT")) + 
  theme_bw() + labs(x="x",y=expression(Omega[y]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p1
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wx.",wch), col=paste0("dgamma.",wch), shape="fsigTILT")) + 
  theme_bw() + labs(x="x",y=expression(Omega[x]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p2
plot_grid(p1,p2,nrow=2)

## HISTOGRAMS
ggplot(dr) + geom_histogram(aes(Wx.CW, fill=dgamma.CW),binwidth=1e-4) + facet_grid(fsigTILT~.)

##CROSS PLOT
filter(DR, abs(dWy)>0, abs(dWx)>0)%>% 
  ggplot(aes(abs(dWy),abs(dWx), col=fsigTILT))+
  geom_point() + theme_bw() + theme(legend.position="top") +
  labs(x=expression(abs(Delta~Omega[y])), y=expression(abs(Sigma~Omega[x])), col=expression(sigma[tilt])) +
  geom_smooth(method="lm", se=FALSE) +
  scale_y_log10() + scale_x_log10() -> p3

ddply(DR, .(fsigTILT), function(df){
  lm(abs(dWx)~abs(dWy), df, abs(dWy)>0) -> m
  coef(summary(m))[1:2, 1:2]-> mt
  c(sigTILT=df$sigTILT[1], Itct=mt[1,1],SEItct=mt[1,2],Slp=mt[2,1],SESlp=mt[2,2])
}) -> .stats

ggplot(.stats, aes(sigTILT, Itct)) + geom_point() +
  geom_errorbar(aes(ymin=Itct-SEItct, ymax=Itct+SEItct)) +
  theme_bw() + labs(
    x=expression(sigma[tilt]), 
    y=expression(abs(Sigma~Omega[x])~"("~abs(Delta~Omega[y])~"=0 )")
  ) +
  scale_y_log10() + scale_x_log10() -> p4
plot_grid(p3,p4,nrow = 2)
