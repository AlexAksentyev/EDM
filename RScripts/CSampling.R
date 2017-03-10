library(methods)
library(dplyr); library(plyr)
library(data.table)

CSampling = setClass(
  "CSampling",
  slots = c(Type = "character", Freq = "numeric", rerror = "numeric"),
  prototype = list(Type = NULL, Freq=5000, rerror = 3e-2),
  contains = "VIRTUAL"
)
CuSampling = setClass("CuSampling", contains = "CSampling", prototype = list(Type="Uniform"))
CmSampling = setClass(
  "CmSampling", contains = "CSampling",
  slots = c(CMPT="numeric", sglFreqGuess = "numeric"),
  prototype = list(Type="Modulated", CMPT=.42, sglFreqGuess = rnorm(1, 3, .01))
)

setGeneric("simSample", def=function(sampling, signal, time, rerror=sampling@rerror, grow=FALSE) standardGeneric("simSample"))
setGeneric("setValue",def=function(object, value) standardGeneric("setValue"))
setGeneric("smplPts", def=function(sampling, time, ...) standardGeneric("smplPts"))

setMethod(
  f="setValue", signature="CSampling",
  definition=function(object, value){
    for(n in names(value)) slot(object, n) <- value[n]
    object
  }
)
setMethod(
  f="simSample", signature = "CuSampling",
  definition=function(sampling, signal, time, rerror = sampling@rerror, grow=FALSE){
    
    if(length(time) < 2) time <- c(0, time)
    
    aerror <- function(x) rerror * ifelse(grow, exp(-.5*signal@beamLam* x), 1)* signal@Num0*signal@Pol
    
    t1 = seq(time[1], time[2], by = 1/sampling@Freq) #uniform sampling
    
    data.table("Time" = t1, "XSgl" = expectation(signal, t1), "FIDrvt" = fiDer(signal, t1))[,Sgl := XSgl + rnorm(length(t1), sd=aerror(t1))]
  }
)

setMethod(
  f="simSample", signature = "CmSampling",
  definition=function(sampling, signal, time, rerror = sampling@rerror, grow=FALSE){
    
    if(sampling@CMPT > 1) sampling@CMPT <- 1
    if(length(time) < 2) time <- c(0, time)
    
    phi = signal@Phase; w0 = signal@wFreq; lam.decoh = signal@decohLam
    P = signal@Pol; N0 = signal@Num0
    fs = sampling@Freq; wg = sampling@sglFreqGuess
    Tpg = pi/wg
    Dt = sampling@CMPT*pi/signal@wFreq
    
    aerror <- function(x) rerror * ifelse(grow, exp(-.5*signal@beamLam* x), 1)* N0*P
    
    .dum <- function(Time) floor((wg*Time+phi)/2/pi)
    Nstt = .dum(time[1])
    Ntot = .dum(time[2])
    # Nnd = floor((duration*wg + phi)/pi); cat(paste("periods", Nnd, "\n"))
    
    tnu = (2*pi*Nstt:Ntot-phi)/wg; tnu <- tnu[tnu>=0]
    tnd = tnu+pi/wg
    t1 = c(tnu, tnd); t1 <- t1[order(t1)]
    # t1 = Tpg*0:Nnd
    
    t2 = seq(-.5*Dt,.5*Dt, 1/fs)
    t3 = rep(t2,length(t1))+rep(t1,each=length(t2))
    
    data.table("Node" = rep(t1,each=length(t2)),
               "Time" = t3, 
               "XSgl" = expectation(signal, t3), 
               "FIDrvt" = fiDer(signal, t3))[,Sgl := XSgl + rnorm(length(t3), sd=aerror(t3))][Time >=time[1] & Time <= time[2],] %>% 
      setattr("CMPT", sampling@CMPT)
  }
)
setMethod(
  f="smplPts", signature = "CuSampling",
  definition=function(sampling, time, ...){
    
    if(length(time) < 2) time <- c(0, time)
    
    t1 = seq(time[1], time[2], by = 1/sampling@Freq)
    
    return(t1)
  }
)
setMethod(
  f="smplPts", signature = "CmSampling",
  definition=function(sampling, time, ...){
    if(sampling@CMPT > 1) sampling@CMPT <- 1
    if(length(time) < 2) time <- c(0, time)
    
    supplied <- list(...)
    if("Phi"%in%names(supplied)) phi = supplied$Phi
    fs = sampling@Freq; wg = sampling@sglFreqGuess
    Dt = sampling@CMPT*pi/wg
    
    .dum <- function(Time) floor((wg*Time+phi)/2/pi)
    Nstt = .dum(time[1])
    Ntot = .dum(time[2])
  
    tnu = (2*pi*Nstt:Ntot-phi)/wg; tnu <- tnu[tnu>=0]
    tnd = tnu+pi/wg
    t1 = c(tnu, tnd); t1 <- t1[order(t1)]
    
    t2 = seq(-.5*Dt,.5*Dt, 1/fs)
    t3 = rep(t2,length(t1))+rep(t1,each=length(t2))
    return(list("Node"=t1, "CMPT"=Dt, "Pts"=t3))
  }
)