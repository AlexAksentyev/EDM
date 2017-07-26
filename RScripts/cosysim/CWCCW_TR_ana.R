library(readr)
library(data.table)
library(plyr); library(dplyr)
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

.degroup <- function(scanvec){
  rle(!is.na(scanvec))->x
  gn=length(x$lengths[!x$values])
  gl <- x$lengths[x$values][1]
  f <- rep(1:gn, each=gl)
  scanvec<-scanvec[!is.na(scanvec)]
  
  split(scanvec,f)
}

.read_data <- function(filename, directory = "~/git/COSYINF/test/"){
  require(reshape2)
  cnt <- switch(filename%%100, "1" = "Sx", "2" = "Sy", "3" = "Sz",
                  "4"="x","5"="a","6"="y","7"="b","8"="l","9"="d")
    
  filename <- paste0(directory,"fort.",filename)
  SaTR = as.numeric(scan(filename,skip=1,na.strings = "g"))
  
  .degroup(SaTR) -> SaTR.list
  ldply(SaTR.list, function(e) data.frame(Turn=e[1], t(e[-(1:2)]))) %>%data.table() -> SaTR
  rm(SaTR.list); SaTR[,.id:=NULL]

  melt.data.table(SaTR, id.vars = "Turn", variable.name = "PID", value.name = cnt) -> SaTR
  SaTR[,Sec:=Turn*1e-6]
  setkey(SaTR,Turn, Sec, PID); SaTR
}

ddt <- function(dt1, dt2, id.vars=c("Turn","PID","Proc","variable")){
  which(id.vars%in%names(dt1)) -> i
  i <- setdiff(1:length(names(dt1)), i)
  id.vars <- id.vars[id.vars!="Proc"]
  setkeyv(dt1, id.vars)
  setkeyv(dt2, id.vars)
  
  merge(dt1,dt2)[,Diff:=value.x-value.y][]
}
## Comparing the horizontally rotated vs unrotated TR  ####
if(FALSE){

  PID <- read_table("~/git/COSYINF/test/fort.100", col_names=c("x","y","d"))
  PID$PID <- paste0("X",1:nrow(PID))
  
  llply(101:109, .read_data) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL
  DL[,Proc:="TR"]
  llply(201:209, .read_data) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL1
  DL1[,Proc:="TR1"]
  iv = c("Turn","Sec","PID","Proc")
  DDL <- ddt(melt(DL,id.vars = iv), melt(DL1,id.vars = iv));
  DL <- rbind(DL,DL1)%>%melt(id.vars=iv)
  
  ggplot(DL, aes(Sec, value, col=Proc)) + 
    geom_point(size=.1) + 
    facet_wrap(~variable,scales="free_y") + 
    scale_color_manual(values=c("red","blue"), breaks=c("TR1","TR"), labels=c("present","absent"), name="XZ-rotation") +
    theme(legend.position="top")

  ggplot(DL[PID%in%c("X1","X2","X3","X4")&Turn<500*300], aes(Turn, l, col=Proc)) + 
    geom_point() + geom_line() +
    theme_bw() + facet_grid(PID~., scales = "free_y") 
}

## Analysis tilt ####

if(FALSE){
  PID <- read_table("~/git/COSYINF/test/fort.100", col_names=c("x","y","d"))
  PID$PID <- paste0("X",1:nrow(PID))
  
  ldply(1:5, function(p) {
    l1 <- 100 + (p-1)*100 + 1:9
    l2 <- 600 + (p-1)*100 + 1:9
    
    llply(l1, .read_data) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL
    DL[,`:=`(Proc="TR",TILT=p*1e-4)]
    llply(l2, .read_data) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL1
    DL1[,`:=`(Proc="TR1",TILT=p*1e-4)]
    iv = c("Turn","Sec","PID","Proc","TILT")
    DDL <- ddt(melt(DL,id.vars = iv), melt(DL1,id.vars = iv));
    DL <- rbind(DL,DL1)%>%melt(id.vars=iv)
  }) %>% data.table() -> DL
  
  ggplot(DL[variable%in%c("Sy","x","a","y","b")],aes(Sec, value, col=as.factor(TILT))) +
    geom_point(size=.2) + theme_bw() +
    facet_grid(variable~Proc,scales="free_y") + 
    # scale_color_manual(values=c("red","blue"), breaks=c("TR1","TR"), labels=c("present","absent"), name="XZ-rotation") +
    theme(legend.position="top")
  
  ggplot(DL, aes(Sec, value, col=Proc)) + 
    geom_point(size=.1) + 
    facet_wrap(~variable,scales="free_y") + 
    scale_color_manual(values=c("red","blue"), breaks=c("TR1","TR"), labels=c("present","absent"), name="XZ-rotation") +
    theme(legend.position="top")
  
  
}

## Analysis solenoid By ####
if(FALSE){
  .loadData(10,11) -> DR
  DR[,`:=`(fBsol=as.factor(Bsol/10),
           # fd=factor(d,labels=formatC(unique(d),0,format="e"))
           fd = cut(d,9,labels=1:9)
  )
  ]
  ggplot(DR[x==1e-8],aes(Wy,Wx,col=fd))+geom_point() + 
    facet_grid(fBsol~Drn,scales="free") +
    theme_bw() #+ theme(legend.position="top")
  
}

## testing Rotation ####
Rotate <- function(invec, angle){
  Ry = c(cos(angle), 0, -sin(angle),
         0, 1, 0,
         sin(angle), 0 , cos(angle))
  
  outvec <- array(3)
  for(i in 1:3) outvec[i]<- sum(Ry[3*(i-1) + 1:3]*invec)
  
  outvec
}

ang <- scan("~/git/COSYINF/test/fort.204",skip=1)
plot(ang[seq(1,length(ang),length.out = 300)], type="l")
pacf(ang)
ll <- length(ang)
s <- array(dim=c(ll,3), dimnames = list(NULL, c("Sx","Sy","Sz"))); s[1,] <- c(0,0,1)
for(n in 2:ll) s[n,] <- Rotate(s[n-1,],ang[n-1])
plot(s[,"Sx"],s[,"Sz"])
plot(s[,"Sz"])

## Comparing the spin maps from TR and TR1 ####

smTR <- scan("~/git/COSYINF/test/fort.111",skip=1)
id = 1; npid = nrow(PID)
smTR <- matrix(smTR[seq(id+1,length(smTR), by=npid+1)], nrow=3, ncol=3, byrow=TRUE)
smTR1 <- scan("~/git/COSYINF/test/fort.212",skip=1,na.strings = "g")
smTR1 <- matrix(smTR1[seq(id+1,length(smTR1), by=npid+1)], nrow=3, ncol=3, byrow=TRUE)

smTR == smTR1

Sn = array(dim=3)
for(i in 1:3) Sn[i] <- sum(smTR[i,]*c(0,0,1))
Sn

### Spectral analysis ####
par(mfrow=c(5,2))
laply(1:5, function(p){
  DL[Proc=="TR1"&variable=="Sy"&TILT==p*1e-4]->sdat
  tsSY <- ts(sdat$value, start=0, end=sdat[nrow(sdat),Sec], deltat=as.numeric(sdat$Sec[2]-sdat$Sec[1]))
  plot(tsSY, main="Tilt, rotated TR", ylab="Sy")
  spec.pgram(tsSY,log="no") -> sps
  sps$freq[which.max(sps$spec)]
})

