library(methods)

CSignal = setClass(
  "CSignal",
  slots = c(
    wFreq0="numeric", wFreq="numeric", 
    Phi0="numeric", Phi="numeric",
    G="numeric"
  ),
  prototype = list(wFreq0=3, Phi0=0, G=7e2)
)

setGeneric("Phase", def=function(object, at) standardGeneric("Phase"))
setMethod(
  f="Phase", signature = "CSignal",
  definition = function(object, at){
    object@wFreq%o%at + object@Phi
  }
)
setGeneric("Signal", def=function(object, at) standardGeneric("Signal"))
setMethod(
  f="Signal", signature="CSignal",
  definition=function(object,at){
    return(data.frame(
      Time=at, Sgl=colSums( sin(Phase(object, at)) )
    ))
  }
)
setGeneric("popPS", def=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2) standardGeneric("popPS"))
setMethod(
  f="popPS", signature = "CSignal",
  definition=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2){
    slot(object, "wFreq") <- object@wFreq0 + object@G*rnorm(Npart, sd=SDdy)^2
    slot(object, "Phi") <- rnorm(Npart, object@Phi0, SDphi)
    
    return(object)
  }
)