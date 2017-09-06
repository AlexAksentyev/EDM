library(readr)
library(data.table)
library(plyr); library(dplyr)
library(ggplot2)
library(cowplot)
library(lattice)

rm(list=ls(all=TRUE))

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


## Analysis ####

if(TRUE){
  
  ldply(1:9, function(p) {
    l2 <- 600 + (p-1)*100 + 1:9
    scan(paste0(from,"fort.1"),skip = 1) -> p1list
    scan(paste0(from,"fort.1"),what=character(), nlines=1)%>%paste0(collapse="_") -> p1name
    scan(paste0(from,"fort.2"),skip = 1)-> p2list
    scan(paste0(from,"fort.2"),what=character(), nlines=1)%>%paste0(collapse="_") -> p2name
    p1 <- p1list[p]
    p2 <- p2list[p] 
    
    llply(l2, read_data, directory=from) %>% Reduce(function(dt1, dt2) merge(dt1,dt2),.) -> DL1
    DL1[,Rotated:="Yes"][,(p1name):=p1][,(p2name):=p2]
    iv = c("Turn","Sec","PID","Rotated",p1name,p2name)
    melt(DL1,id.vars = iv)
  }) %>% data.table() -> DL
  
  DL[variable=="Sy"&Turn==2] -> sdat
  sdat[,.(DeltaSy=value,DT=magnet_length/(.4855*2.9979e3)),by=c("magnet_length","magnet_strength","PID")] ->df
  data.table(merge(PID,df))->df
  ggplot(df[PID%in%paste0("X",c(1:9))], aes(magnet_length,DeltaSy, col=as.factor(x),shape=as.factor(y))) + 
    geom_point(size=2) + geom_line()+ #scale_y_log10() + 
    theme_bw() + #theme(legend.position="top") +
    facet_grid(magnet_strength~.,scales = "free_y") +
    scale_x_continuous(labels=function(x)formatC(x,2,format="e"))
  
  
  ggplot(df,aes(magnet_length, DeltaSy, col=as.factor(magnet_strength))) + 
    geom_point() + geom_line() + 
    theme_bw() + theme(legend.position="top",axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_grid(y~x,scales = "free_y") + 
    scale_x_continuous(labels=function(x)formatC(x,2,format="e"))
  
  
  ggplot(DL[variable%in%c("Sy","Sx","x","y")&Turn<100&PID%in%paste0("X",c(3,8))&magnet_strength==.46],
         aes(Turn, value, col=PID)) +
    geom_line(size=.3) + geom_point(size=1) + theme_bw() + #scale_y_log10()+
    facet_grid(variable~magnet_length,scales="free_y") + 
    theme(legend.position="top") #+ scale_color_manual(values=c("blue","red"))
  
  ## spectral analysis
  if(FALSE){
    par(mfrow=c(3,2))
    laply(1:3, function(p){
      DL[Rotated=="Yes"&PID=="X1"&variable=="Sy"&TILT==unique(TILT)[p]]->sdat
      tsSY <- ts(sdat$value, start=0, end=sdat[nrow(sdat),Sec], deltat=as.numeric(sdat$Sec[2]-sdat$Sec[1]))
      plot(tsSY, ylab="Sy", main=sdat$TILT[1])
      spec.pgram(tsSY,log="no") -> sps
      sps$freq[which.max(sps$spec)]*2*pi
    })
    par(mfrow=c(1,1))
  }

  print(PID)
  
}


