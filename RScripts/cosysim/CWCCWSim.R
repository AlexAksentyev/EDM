library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(lattice)
library(scatterplot3d)
library(mosaic)

.complete <- function(DT){
  Nrev = attr(DT,"Parameters")$Nrev #revolutions per simulation
  Nrps = 1e6 #revolutions per second
  
  DT[,':='(
    ThXY = atan(Sx/Sy),
    ThXZ = atan(Sx/Sz),
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
  # nbar = read_csv(ddf, skip=2, n_max=1); names(nbar) <- c("x","y","z")
  Direct <- data.table(read_csv(ddf, skip=4))
  attr(Direct, "Parameters")<-para
  # attr(Direct, "Nbar")<-nbar
  para =  read_csv(rdf, n_max=1)
  # nbar = read_csv(rdf, skip=2, n_max=1); names(nbar) <- c("x","y","z")
  Reverse <- data.table(read_csv(rdf, skip=4))
  attr(Reverse, "Parameters")<-para
  # attr(Reverse, "Nbar")<-nbar
  
  DR <- .prep(Direct, Reverse)
  
  st0 <- formatC(para$muTILT, 2, format = "e")
  st1 <- formatC(para$sigTILT, 2, format = "e")
    
  DR[,`:=`(muTILT=st0, sigTILT=st1, GSDP=para$GSDP, GSFP=para$GSFP)]
  DR
}) %>% data.table() -> DR

ldply(1:6, function(k){
  ddf = paste0("~/git/COSYINF/test/fort.",(6+2*k))
  rdf = paste0("~/git/COSYINF/test/fort.",(6+2*k+1))
  
  nbarD = read_csv(ddf, skip=2, n_max=1); names(nbarD) <- c("x","y","z")
  nbarR = read_csv(rdf, skip=2, n_max=1); names(nbarR) <- c("x","y","z")
  rbind(mutate(nbarD, Dir="CW"), mutate(nbarR, Dir="CCW"))
}) -> Nbardf

##X-DGAMMA OMEGA PLOT 
dr <- DR[,.(sigTILT, sigTILT, dgamma.CW, dgamma.CCW, x.CW, x.CCW, Wy.CW, Wy.CCW, Wx.CW, Wx.CCW)]
wch = "CCW"
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wy.",wch), col=paste0("dgamma.",wch), shape="sigTILT")) + 
  theme_bw() + scale_color_continuous(high = "red", low="blue") + 
  labs(x="x",y=expression(Omega[y]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p1
ggplot(dr) + geom_point(aes_string(paste0("x.",wch), paste0("Wx.",wch), col=paste0("dgamma.",wch), shape="sigTILT")) + 
  theme_bw() + scale_color_continuous(high = "red", low="blue") + 
  labs(x="x",y=expression(Omega[x]), col=expression(Delta~gamma),shape=expression(sigma[tilt])) -> p2
plot_grid(p1,p2,nrow=2)

## HISTOGRAMS
ggplot(dr) + geom_histogram(aes(Wx.CW, fill=dgamma.CW),binwidth=1e-4) + facet_grid(sigTILT~.)

##CROSS PLOT
filter(DR, abs(dWy)>0, abs(dWx)>0)%>% 
  ggplot(aes(abs(dWy),abs(dWx), col=sigTILT))+
  geom_point() + theme_bw() + theme(legend.position="top") +
  labs(x=expression(abs(Delta~Omega[y])), y=expression(abs(Sigma~Omega[x])), col=expression(sigma[tilt])) +
  geom_smooth(method="lm", se=FALSE) +
  scale_y_log10() + scale_x_log10() -> p3

ddply(DR, .(sigTILT), function(df){
  lm(abs(dWx)~abs(dWy), df, abs(dWy)>0) -> m
  coef(summary(m))[1:2, 1:2]-> mt
  c(Itct=mt[1,1],SEItct=mt[1,2],Slp=mt[2,1],SESlp=mt[2,2])
}) -> .stats

ggplot(.stats, aes(as.numeric(sigTILT), Itct)) + geom_point() +
  geom_errorbar(aes(ymin=Itct-SEItct, ymax=Itct+SEItct)) +
  theme_bw() + labs(
    x=expression(sigma[tilt]), 
    y=expression(abs(Sigma~Omega[x])~"("~abs(Delta~Omega[y])~"=0 )")
  ) +
  scale_y_log10() + scale_x_log10() -> p4
plot_grid(p3,p4,nrow = 2)

## injection point analysis
dr <- DR[x.CW==x.CCW&dgamma.CW==dgamma.CCW,][,slp:=dWx/dWy][,.(x.CW, dgamma.CW, dWx, dWy, dWz, slp, sigTILT)]
dr[,sslp:=as.vector(scale(slp,center=.5*(max(slp)+min(slp)),scale=ifelse(max(slp)-min(slp)!=0,max(slp)-min(slp),1))),by=sigTILT]
dr[sigTILT=="1.00e-04",] %>% 
  ggplot(aes(x.CW, dgamma.CW)) +
  geom_contour(aes(z=sslp,col=..level..)) +
  scale_color_continuous(low="blue",high="red") + theme_bw()

levelplot(sslp~x.CW*dgamma.CW, dr[sigTILT=="1.00e-04",])

ggplot(dr[x.CW==-.001&sigTILT!="0.00e+00"], aes(dgamma.CW,sslp)) + geom_point(aes(col=sigTILT))
ggplot(dr[dgamma.CW==1e-5&sigTILT=="1.00e-06"], aes(x.CW,sslp)) + geom_point(aes(col=sigTILT))

dr[,
   pcol:=derivedFactor(
     "red" = sslp>.90*max(sslp), 
     "magenta" = sslp>.60*max(sslp),
     "green" = sslp>.20*max(sslp), 
     .default="blue", 
     .method="first"
   )]
with(dr[sigTILT=="1.00e-05"],{scatterplot3d(
  x.CW, dgamma.CW, sslp,
  color=pcol, pch=19, type="o",
  highlight.3d = TRUE
)})


## ONE DIRECTION ####

DR[x.CW==x.CCW&dgamma.CW==dgamma.CCW,
   .(x.CW, dgamma.CW, sigTILT, GSDP, GSFP, Wx.CW, Wx.CCW)
]->dr

ggplot(dr[sigTILT=="1.00e-03"], aes(x.CW, dgamma.CW)) + geom_tile(aes(fill=Wx.CW))+
  scale_fill_continuous(high="red")
with(dr[sigTILT=="1.00e-03"],{scatterplot3d(
  x.CW, dgamma.CW, Wx.CW, 
  type="o", highlight.3d = TRUE,
  xlab="x",ylab=expression(Delta~gamma)
)})

ggplot(dr[sigTILT=="1.00e-03"&x.CW==.001], aes(dgamma.CW,Wx.CW)) + geom_point() -> p1
ggplot(dr[sigTILT=="1.00e-03"&dgamma.CW==1e-5], aes(x.CW,Wx.CW)) + geom_point() -> p2
plot_grid(p1,p2,nrow=2)

lm(Wx.CW~dgamma.CW, dr[sigTILT=="1.00e-03"&x.CW==.001]) %>% summary %>% coef
lm(Wx.CW~x.CW, dr[sigTILT=="1.00e-03"&dgamma.CW==1e-5&x.CW<=0]) %>% summary %>% coef

ddply(dr, "sigTILT", function(df){
  ddply(df, "x.CW", function(sub){
   coef(summary(lm(Wx.CW~dgamma.CW, sub))) -> mod
   c("Estimate" = mod[2,1], "SE" = mod[2,2])
  }) %>% melt(id.vars=c("Estimate","SE"), variable.name="Cntrl", value.name="At")-> df1
  ddply(df, "dgamma.CW", function(sub){
    coef(summary(lm(Wx.CW~x.CW, sub[sub$x.CW<=0,]))) -> mod
    c("Estimate" = mod[2,1], "SE" = mod[2,2])
  }) %>% melt(id.vars=c("Estimate","SE"), variable.name="Cntrl", value.name="At") -> df2
  df1<-mutate(df1, GSDP = df$GSDP[1], GSFP = df$GSFP[1])
  df2<-mutate(df2, GSDP = df$GSDP[1], GSFP = df$GSFP[1])
  rbind(df1,df2)
}) -> slpdf

dlply(slpdf, "sigTILT",function(df){
  ggplot(df) + 
    geom_histogram(aes(Estimate, fill=Cntrl), show.legend = FALSE) + 
    geom_vline(xintercept = 0) +
    theme_bw() + ggtitle(bquote(sigma[tilt]~"="~.(df$sigTILT[1])))
}) -> pl
plot_grid(plotlist=pl,nrow=3,ncol=2)

DR[,
  .(CW = mean(Wx.CW), CCW = mean(Wx.CCW), SCW = sd(Wx.CW), SCCW = sd(Wx.CCW)),
  by=muTILT
] -> x0
laply(unique(DR$muTILT)%>%as.numeric, function(mu){th=rnorm(1000,mu,1e-4);s=sin(th);mean(s)}) -> m
x0$Model <- m
x0[,`:=`(CW = CW/CW[2], CCW=CCW/CCW[2], Model=Model/Model[2])]
x0
