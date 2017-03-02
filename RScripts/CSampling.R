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
  prototype = list(Type="Modulated", CMPT=.33, sglFreqGuess = rnorm(1, 3, .01))
)

setGeneric("simSample", def=function(sampling, signal, time, rerror=sampling@rerror) standardGeneric("simSample"))
setGeneric("setValue",def=function(object, value) standardGeneric("setValue"))

setMethod(
  f="setValue", signature="CSampling",
  definition=function(object, value){
    for(n in names(value)) slot(object, n) <- value[n]
    object
  }
)
setMethod(
  f="simSample", signature = "CuSampling",
  definition=function(sampling, signal, time, rerror = sampling@rerror){
    
    aerror <- rerror * signal@Num0*signal@Pol
    
    t1 = seq(time[1], time[2], by = 1/sampling@Freq) #uniform sampling
    
    data.table("Time" = t1, "XSgl" = expectation(signal, t1), "FIDrvt" = fiDer(signal, t1))[,Sgl := XSgl + rnorm(length(t1), sd=aerror)]
  }
)

setMethod(
  f="simSample", signature = "CmSampling",
  definition=function(sampling, signal, time, rerror = sampling@rerror){
    
    if(sampling@CMPT > 1) sampling@CMPT <- 1
    
    phi = signal@Phase; w0 = signal@wFreq; lam.decoh = signal@decohLam
    P = signal@Pol; N0 = signal@Num0
    fs = sampling@Freq; wg = sampling@sglFreqGuess
    Tpg = pi/wg
    Dt = sampling@CMPT*pi/signal@wFreq
    
    aerror <- rerror * N0*P
    
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
               "FIDrvt" = fiDer(signal, t3))[,Sgl := XSgl + rnorm(length(t3), sd=aerror)][Time >=time[1] & Time <= time[2],] %>% 
      setattr("CMPT", sampling@CMPT)
  }
)
