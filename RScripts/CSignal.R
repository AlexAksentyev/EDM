library(methods)
library(ggplot2)

CSignal = setClass(
  "CSignal",
  slots = c(wFreq0="numeric", Phi0="numeric",
            G="numeric", SDdy="numeric", SDphi="numeric",
            PS="data.frame",
            Pol="data.frame"),
  prototype = list(wFreq0=3, Phi0=0, G=7e2)
)

setGeneric("Phase", def=function(object, at) standardGeneric("Phase"))
setMethod(
  f="Phase", signature = "CSignal",
  definition = function(object, at){
    object@PS$wFreq%o%at + object@PS$Phi
  }
)
setGeneric("Signal", def=function(object, at) standardGeneric("Signal"))
setMethod(
  f="Signal", signature="CSignal",
  definition=function(object,at){
    object@Pol <- data.frame(
      Time=at, Sgl=colSums( sin(Phase(object, at)) )
    )
    return(object)
  }
)
setGeneric("popPS", def=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2) standardGeneric("popPS"))
setMethod(
  f="popPS", signature = "CSignal",
  definition=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2){
    object@SDdy <- SDdy
    object@SDphi <- SDphi
    
    df = data.frame(
      wFreq = object@wFreq0 + object@G*rnorm(Npart, sd=SDdy)^2,
      Phi = rnorm(Npart, object@Phi0, SDphi)
    )
    attr(df$wFreq, "Synch") <- object@wFreq0
    attr(df$wFreq, "SD") <- object@G*object@SDdy^2
    attr(df$Phi, "Synch") <- object@Phi0
    attr(df$Phi, "SD") <- object@SDphi
    
    slot(object, "PS") <- df
    
    return(object)
  }
)
