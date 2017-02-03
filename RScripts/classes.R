library(methods)
library(dplyr); library(plyr)

#### structural model ####

CModel = setClass(
  "CModel",
  slots = c(wFreq="numeric", Phase="numeric", Num0="numeric", beamLam = "numeric", Pol="numeric", decohLam="numeric"),
  prototype = list(wFreq=3, Phase=0, Num0=6730, beamLam=-1/2000, Pol=.4, decohLam=log(.25)/1000)
)

setGeneric("setValue",def=function(object, value) standardGeneric("setValue"))
setGeneric("expectation", def=function(object, at) standardGeneric("expectation"))
setGeneric("fiDer", def=function(object, at) standardGeneric("fiDer"))
setGeneric("timeDer", def=function(object, at) standardGeneric("timeDer"))

#### Sampling types ####
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

setGeneric("simSample", def=function(sampling, signal, duration, rerror=sampling@rerror) standardGeneric("simSample"))

#### method definitions ####
setMethod(
  f="setValue", signature="CModel", 
  definition=function(object, value){
    for(n in names(value)) slot(object, n) <- value[n]
    object
  }
)
setMethod(
  f="expectation", signature="CModel",
  definition=function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam; lam.beam <- object@beamLam
    
    N0 * exp(lam.beam*at) * (1 + P*exp(lam.decoh*at)*sin(w*at + phi))
  }
)
setMethod(
  f="fiDer", signature="CModel",
  definition=function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam; lam.beam <- object@beamLam
    
    N0*exp(lam.beam*at) * P*exp(lam.decoh*at) * cos(w*at + phi)
  }
)
setMethod(
  f="timeDer", signature = "CModel",
  definition = function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam; lam.beam <- object@beamLam
    
    etd = exp(lam.decoh*at)
    etb = exp(lam.beam*at)
    
    lam.beam*N0*etb * (1 + P*etd*sin(w*at + phi)) +
      N0*etb * P*(lam.decoh*etd*sin(w*at + phi) + w*etd*cos(w*at + phi))
  }
)

setMethod(
  f="setValue", signature="CSampling",
  definition=function(object, value){
    for(n in names(value)) slot(object, n) <- value[n]
    object
  }
)
setMethod(
  f="simSample", signature = "CuSampling",
  definition=function(sampling, signal, duration, rerror = sampling@rerror){
    
    aerror <- rerror * signal@Num0*signal@Pol
    
    t1 = seq(0, duration, by = 1/sampling@Freq) #uniform sampling
    
    data.frame("Time" = t1, "XSgl" = expectation(signal, t1), "FIDrvt" = fiDer(signal, t1)) %>%  
      mutate(Sgl = XSgl + rnorm(length(t1), sd=aerror))
  }
)

setMethod(
  f="simSample", signature = "CmSampling",
  definition=function(sampling, signal, duration, rerror = sampling@rerror){
    
    phi = signal@Phase; w0 = signal@wFreq; lam.decoh = signal@decohLam
    P = signal@Pol; N0 = signal@Num0
    fs = sampling@Freq; wg = sampling@sglFreqGuess
    Tpg = pi/wg
    Dt = sampling@CMPT
    
    aerror <- rerror * N0*P
    
    Nnd = floor((duration*wg + phi)/pi); cat(paste("periods", Nnd, "\n"))
    t1 = Tpg*0:Nnd
    t2 = seq(-.5*Dt,.5*Dt, 1/fs)
    t3 = rep(t2,length(t1))+rep(t1,each=length(t2))
    
    data.frame("Node" = rep(t1,each=length(t2)),"Time" = t3, "XSgl" = expectation(signal, t3), "FIDrvt" = fiDer(signal, t3)) %>% 
      `attr<-`("CMPT", sampling@CMPT) %>% 
      mutate(Sgl = XSgl + rnorm(length(t3), sd=aerror)) %>%
      filter(Time >=0 & Time <= duration)
  }
)
