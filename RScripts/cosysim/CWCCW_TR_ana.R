library(readr)
library(data.table)
library(plyr); library(dplyr)
library(ggplot2)
library(cowplot)
library(lattice)

## DEFINITIONS ####

degroup <- function(scanvec){
  rle(!is.na(scanvec))->x
  gn=length(x$lengths[!x$values])
  gl <- x$lengths[x$values][1]
  f <- rep(1:gn, each=gl)
  scanvec<-scanvec[!is.na(scanvec)]
  
  split(scanvec,f)
}

read_data <- function(filename, directory = "~/git/COSYINF/BNL/"){
  require(reshape2)
  cnt <- switch(filename%%100, "1" = "Sx", "2" = "Sy", "3" = "Sz",
                  "4"="x","5"="a","6"="y","7"="b","8"="l","9"="d")
    
  filename <- paste0(directory,"fort.",filename)
  SaTR = as.numeric(scan(filename,skip=1,na.strings = "g"))
  
  degroup(SaTR) -> SaTR.list
  ldply(SaTR.list, function(e) data.frame(Turn=e[1], t(e[-(1:2)]))) %>%data.table() -> SaTR
  rm(SaTR.list); SaTR[,.id:=NULL]

  melt.data.table(SaTR, id.vars = "Turn", variable.name = "PID", value.name = cnt) -> SaTR
  SaTR[,Sec:=Turn*1e-6]
  setkey(SaTR,Turn, Sec, PID); SaTR
}

ddt <- function(dt1, dt2, id.vars=c("Turn","PID","Rotated","variable")){
  which(id.vars%in%names(dt1)) -> i
  i <- setdiff(1:length(names(dt1)), i)
  id.vars <- id.vars[id.vars!="Rotated"]
  setkeyv(dt1, id.vars)
  setkeyv(dt2, id.vars)
  
  merge(dt1,dt2)[,Diff:=value.x-value.y][]
}

## ANALYSIS ####
from = "~/git/COSYINF/test/"

## reading particle initial data ####
PID <- read_table(paste0(from,"fort.100"), col_names=c("x","y","d"))
PID$PID <- paste0("X",1:nrow(PID))

## Analysis tilt ####

if(FALSE){
  
  ldply(1:2, function(p) {
    l1 <- 100 + (p-1)*100 + 1:9
    l2 <- 600 + (p-1)*100 + 1:9
    tiltp <- (p-1)*0.0057
    
    llply(l1, read_data, directory=from) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL
    DL[,`:=`(Rotated="No",TILT=tiltp)]
    llply(l2, read_data, directory=from) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL1
    DL1[,`:=`(Rotated="Yes",TILT=tiltp)]
    iv = c("Turn","Sec","PID","Rotated","TILT")
    DDL <- ddt(melt(DL,id.vars = iv), melt(DL1,id.vars = iv));
    DL <- rbind(DL,DL1)%>%melt(id.vars=iv)
    # melt(DL,id.vars = iv)
  }) %>% data.table() -> DL
  
  ggplot(DL[variable%in%c("Sy","Sx","Sz","x","y")],aes(Sec, value, col=Rotated)) +
    geom_line(size=.2) + theme_bw() +
    facet_grid(variable~TILT,scales="free_y") + 
    theme(legend.position="top")
  
  sdat[Rotated=="Yes",.(Th = asin(value), Sec, PID), by=c("TILT","variable","PID")][,.(W=(Th-Th[1])/Sec, Sec), by=c("TILT","variable","PID")]->Wdat
  ggplot(Wdat[Sec>.25&variable=="Sy"], aes(Sec, W, col=PID)) + geom_line() + 
    facet_grid(TILT~variable,scales="free_y") + geom_smooth(method="gam") + 
    theme(legend.position="top")
  
  ## spectral analysis
  par(mfrow=c(2,2))
  laply(1:2, function(p){
    DL[Rotated=="Yes"&variable=="Sy"&TILT==(p-1)*0.0057]->sdat
    tsSY <- ts(sdat$value, start=0, end=sdat[nrow(sdat),Sec], deltat=as.numeric(sdat$Sec[2]-sdat$Sec[1]))
    plot(tsSY, ylab="Sy")
    spec.pgram(tsSY,log="no") -> sps
    sps$freq[which.max(sps$spec)]
  })
  
}


