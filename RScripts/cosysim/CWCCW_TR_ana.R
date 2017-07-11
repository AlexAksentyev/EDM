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
  
  llply(filenames, function(k){
    filename = paste0("~/git/COSYINF/test/fort.",k)
    para = read_csv(filename, n_max=1)
    DT <- data.table(read_csv(filename, skip=2))
    setattr(DT, "Parameters", para)
    setkey(DT, muTILT, sigTILT, Trl, x, a, y, b, l, d)
    .complete(DT)
  }) 
  
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
  SyTR = scan(filename, " ",skip=1)
  i = which(SyTR==SyTR[1])
  f = rep(1:length(i), each=diff(i)[1])
  split(SyTR,f) -> SyTR.list
  ldply(SyTR.list, function(e) data.frame(t(as.numeric(e[-1]))), .id="Turn") %>%data.table() -> SYTR
  melt.data.table(SYTR, id.vars = "Turn", variable.name = "Part", value.name = "Sy")
}

## Comparing the horizontally rotated vs unrotated TR  ####
fileTR = "~/git/COSYINF/test/fort.7"
fileTR1 = "~/git/COSYINF/test/fort.8"

SYTR <- .readSpin(fileTR)
SYTR1 <- .readSpin(fileTR1)

SYTR[,Proc:="TR"]; SYTR1[,Proc:="TR1"]
rbind(SYTR,SYTR1) -> SYTRb
SYTRb[,Turn:=as.numeric(Turn)]
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

## Analysis ####

.loadData(10) -> DR

DR2 <- rbind(DR$FWD10[,Dir:="FWD"],DR$REV11[,Dir:="REV"])
ggplot(DR2[y==1e-8&d==min(abs(d))],aes(x,Wx,col=Dir)) + geom_point()+geom_line() +
  theme(legend.position="top") 
