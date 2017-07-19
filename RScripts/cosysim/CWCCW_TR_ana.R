library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(lattice)

## DEFINITIONS ####

.complete <- function(DT){
  Nrev = attr(DT,"Parameters")$Nrev #revolutions per simulation
  Nrps = attr(DT,"Parameters")$Nrps  #revolutions per second
  
  DT[,`:=`(
    ThXY = atan(Sx/Sy),
    ThXZ = atan(Sx/Sz),
    ThYZ = atan(Sy/Sz)
  )
  ]
  DT[,`:=`(
    Wz = ThXY/Nrev*Nrps,
    Wy = ThXZ/Nrev*Nrps,
    Wx = ThYZ/Nrev*Nrps,
    T.Sy = 2*pi/ThYZ*Nrev
  )]
  
}

.prep <- function(DT, RT){
  .complete(DT)
  .complete(RT)
  

  
  # names(DT) <- paste0(names(DT),".FWD")
  # names(RT) <- paste0(names(RT),".REV")
  
  DT[rep(1:nrow(DT),each=nrow(RT)),]->DT
  RT[rep(1:nrow(RT),nrow(RT)),]->RT
  cbind(Direct, RT) -> DR
  DR[,':='(
    dWx = Wx.FWD+Wx.REV,
    dWy = Wy.FWD-Wy.REV,
    dWz = Wz.FWD+Wz.REV
  )]
  
  return(DR)
}

.loadData <- function(...){
  
  filenames = c(...)
  n = floor(length(filenames)/2)
  n <- ifelse(n>0,n,1)
  names(filenames) <- ifelse(filenames%%2,"REV","FWD")%>%paste0(filenames)
  
  ldply(filenames, function(k){
    filename = paste0("~/git/COSYINF/test/fort.",k)
    drn <- ifelse(k%%2,"REV","FWD")%>%paste0(k)
    para = read_csv(filename, n_max=1)
    DT <- data.table(read_csv(filename, skip=2))
    setattr(DT, "Parameters", para)
    # setkey(DT, muTILT, sigTILT, Trl, x, a, y, b, l, d)
    .complete(DT)
  }, .id="Drn") %>%data.table()
  
  # filenameD = paste0("~/git/COSYINF/test/fort.",k)
  # filenameR = paste0("~/git/COSYINF/test/fort.",k+1)
  # 
  # para = read_csv(filenameD, n_max=1)
  # DT <- data.table(read_csv(filenameD, skip=2))
  # setattr(DT, "Parameters", para)
  # setkey(DT, muTILT, sigTILT, Trl, x, a, y, b, l, d)
  # para = read_csv(filenameR, n_max=1)
  # RT <- data.table(read_csv(filenameR, skip=2))
  # setattr(RT, "Parameters", para)
  # setkey(RT, muTILT, sigTILT, Trl, x, a, y, b, l, d)
  
  # .complete(DT)
  # .complete(RT)
  # 
  # list(FWD=DT, REV=RT)
  
}

.readSpin <- function(filename){
  require(reshape2)
  cnt <- switch(filename%%100, "1" = "Sx", "2" = "Sy", "3" = "Sz")
  filename <- paste0("~/git/COSYINF/test/fort.",filename)
  SaTR = as.numeric(scan(filename,skip=1,na.strings = "g"))
  
  rle(!is.na(SaTR))->x
  gn=length(x$lengths[!x$values])
  gl <- x$lengths[x$values][1]
  f <- rep(1:gn, each=gl)
  SaTR<-SaTR[!is.na(SaTR)]
  
  split(SaTR,f) -> SaTR.list
  ldply(SaTR.list, function(e) data.frame(t(e[-1])), .id="Turn") %>%data.table() -> SaTR
  rm(SaTR.list)

  melt.data.table(SaTR, id.vars = "Turn", variable.name = "Part", value.name = cnt) -> SaTR
  setkey(SaTR,Turn, Part); SaTR
}

## Comparing the horizontally rotated vs unrotated TR  ####
SXTR <- .readSpin(101)
SYTR <- .readSpin(102)
SZTR <- .readSpin(103)
SXYZTR <- merge(SXTR,SYTR)%>%merge(SZTR)
SXTR <- .readSpin(201)
SYTR <- .readSpin(202)
SZTR <- .readSpin(203)
SXYZTR1 <- merge(SXTR,SYTR)%>%merge(SZTR)
rm(SXTR,SYTR,SZTR)

SXYZTR[,Proc:="TR"]; SXYZTR1[,Proc:="TR1"]
rbind(SXYZTR,SXYZTR1) -> SXYZTRb
SXYZTRb[,Turn:=as.numeric(Turn)]

ggplot(SXYZTRb,aes(Turn,Sy, col=Proc)) + geom_line() + 
  facet_grid(Proc~.,scales="free_y")
# subset(cbind(SYTR,SYTR1),select=-c(4,5,6,8)) -> SYTRb2
# setnames(SYTRb2, c("Turn","Part","Sy.TR","Sy.TR1"))
# SYTRb2[,DSy:=Sy.TR-Sy.TR1]
# SYTRb2[Turn%in%c("1","50","100"),.(sdDSy=sd(DSy)),by="Turn"]

# ggplot(SYTRb2[Turn%in%c("1","50","100")], aes(DSy)) + geom_histogram() + 
#   facet_grid(.~Turn,scales="free_x") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x="TR - TR1 Sy difference")

ggplot(SYTRb[Part%in%c("X1","X41", "X81")], aes(Turn, Sy, col=Proc)) + 
  geom_line() + geom_point() + 
  theme_bw() + facet_grid(Part~., scales = "free_y")

## Analysis tilt ####

.loadData(10,11) -> DR

DR2 <- rbind(DR$FWD10[,Dir:="FWD"],DR$REV11[,Dir:="REV"])
ggplot(DR2[d==-1e-3,.(meanWx=mean(Wx),SEWx=sd(Wx)/sqrt(10)),by=c("muTILT","x","Dir")],aes(x,meanWx,col=Dir)) + 
  geom_point()+ geom_pointrange(aes(ymin=meanWx-SEWx,ymax=meanWx+SEWx)) +
  geom_line() +
  theme(legend.position="top") + facet_grid(muTILT~.,scales="free_y")

## Analysis solenoid By ####
.loadData(10,11) -> DR
DR[,`:=`(fBsol=as.factor(Bsol/10),
         # fd=factor(d,labels=formatC(unique(d),0,format="e"))
         fd = cut(d,9,labels=1:9)
         )
   ]
ggplot(DR[x==1e-8],aes(Wy,Wx,col=fd))+geom_point() + 
  facet_grid(fBsol~Drn,scales="free") +
  theme_bw() #+ theme(legend.position="top")
