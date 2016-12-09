library(methods)

#### structural model ####

CModel = setClass(
  "CModel",
  slots = c(wFreq="numeric", Phase="numeric",Num0="numeric", Pol="numeric", decohLT="numeric"),
  prototype = list(wFreq=3, Phase=0, Num0=6730, Pol=.4, decohLT=log(.25)/1000)
)

setGeneric("setValue",def=function(object, value) standardGeneric("setValue"))
setGeneric("expectation", def=function(object, at) standardGeneric("expectation"))
setGeneric("fisherInfo", def=function(object, at) standardGeneric("fisherInfo"))

#### Sampling types ####
CSampling = setClass(
  "CSampling",
  slots = c(Type = "character", Freq = "numeric"),
  prototype = list(Type = NULL, Freq=5000),
  contains = "VIRTUAL"
)
CuSampling = setClass("CuSampling", contains = "CSampling", prototype = list(Type="Uniform"))
CmSampling = setClass(
  "CmSampling", contains = "CSampling",
  slots = c(Compaction="numeric", sglFreqGuess = "numeric"),
  prototype = list(Type="Modulated", Compaction=.33, sglFreqGuess = rnorm(1, 3, .001))
)


setGeneric("simSample", def=function(object, signal, duration, rerror=NA) standardGeneric("simSample"))

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
    lam.decoh <- object@decohLT
    
    N0 * (1 + P*exp(lam.decoh*at)*sin(w*at + phi))
  }
)
setMethod(
  f="fisherInfo", signature="CModel",
  definition=function(object, at){
    N0 <- object@Num0; P <- object@Pol; w <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLT
    
    N0*P*exp(lam.decoh*at)*cos(w*at + phi)
  }
)

setMethod(
  f="simSample", signature = "CuSampling",
  definition=function(object, signal, duration, rerror = NA){
    if(is.na(rerror)) rerror <- 3e-2
    
    t1 = seq(0, duration, by = 1/object@Freq) #uniform sampling
    
    data.frame("Time" = t1, "XSgl" = expectation(signal, t1), "fishInfo" = fisherInfo(signal, t1)) %>%  
      mutate(Sgl = XSgl + rnorm(length(t1), sd=rerror*XSgl))
  }
)
setMethod(
  f="simSample", signature = "CmSampling",
  definition=function(object, signal, duration, rerror = NA){
    if(is.na(rerror)) rerror <- 3e-2
    
    phi = signal@Phase; w0 = signal@wFreq; lam.decoh = signal@decohLT
    P = signal@Pol; N0 = signal@Num0
    fs = object@Freq; wg = object@sglFreqGuess
    cptn = object@Compaction
    
    Nprd = round((duration*w0 + phi)/(2*pi)); cat(paste("periods", Nprd, "\n"))
    Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)]
    tn = ((0:(2*Nprd))*pi - phi)/wg; cat(paste("last z-c", tn[length(tn)], "\n"))
    t2 = laply(tn, function(ti) ti+Dt) %>% c 
    
    x = expectation(signal, t2);  DeltaS = P*cptn*exp(lam.decoh*t2) # information condition
    t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0]
    t2 <- t2[order(t2)]
    
    data.frame("Time" = t2, "XSgl" = expectation(signal, t2), "Drvt" = fisherInfo(signal, t2)) %>%
      mutate(Sgl = XSgl + rnorm(length(t2), sd=rerror*XSgl))
  }
)
