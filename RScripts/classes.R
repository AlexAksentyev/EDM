library(methods)
library(dplyr); library(plyr)

#### structural model ####

CModel = setClass(
  "CModel",
  slots = c(wFreq="numeric", Phase="numeric",Num0="numeric", Pol="numeric", decohLam="numeric"),
  prototype = list(wFreq=3, Phase=0, Num0=6730, Pol=.4, decohLam=log(.25)/1000)
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
  slots = c(Compaction="numeric", sglFreqGuess = "numeric"),
  prototype = list(Type="Modulated", Compaction=.33, sglFreqGuess = rnorm(1, 3, .01))
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
    lam.decoh <- object@decohLam
    
    N0 * (1 + P*exp(lam.decoh*at)*sin(w*at + phi))
  }
)
setMethod(
  f="fiDer", signature="CModel",
  definition=function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam
    
    N0*P*exp(lam.decoh*at)*cos(w*at + phi)
  }
)
setMethod(
  f="timeDer", signature = "CModel",
  definition = function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam
    
    et = exp(lam.decoh*at)
    
    N0*P*(lam.decoh*et*sin(w*at + phi) + w*et*cos(w*at + phi))
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
    cptn = sampling@Compaction
    
    aerror <- rerror * N0*P
    
    Nprd = round((duration*w0 + phi)/(2*pi)); cat(paste("periods", Nprd, "\n"))
    Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)]
    tn = ((0:(2*Nprd))*pi - phi)/wg; cat(paste("last z-c", tn[length(tn)], "\n"))
    t2 = laply(tn, function(ti) ti+Dt) %>% c 
    
    x = expectation(signal, t2);  DeltaS = P*cptn*exp(lam.decoh*t2) # information condition
    t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0 & t2 <= duration]
    t2 <- t2[order(t2)]
    
    data.frame("Time" = t2, "XSgl" = expectation(signal, t2), "FIDrvt" = fiDer(signal, t2)) %>% 
      `attr<-`("Compaction", sampling@Compaction) %>% 
      mutate(Sgl = XSgl + rnorm(length(t2), sd=aerror))
  }
)
