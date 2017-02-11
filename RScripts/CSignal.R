library(methods)
library(ggplot2)

CSignal = setClass(
  "CSignal",
  slots = c(
    Synch="numeric", SD="numeric",G="numeric",
    PS="data.frame"
  ),
  prototype = list(Synch=c(wFreq=3, Phi=0), G=7e2)
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
    data.frame(
      Time=at, Sgl=colSums( sin(Phase(object, at)) )
    )
  }
)
setGeneric("popPS", def=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2) standardGeneric("popPS"))
setMethod(
  f="popPS", signature = "CSignal",
  definition=function(object, Npart=1e3, SDdy=1e-3, SDphi=1e-2){
    object@SD <- c(dy = SDdy, phi = SDphi)
    
    df = data.frame(
      wFreq = object@Synch["wFreq"] + object@G*rnorm(Npart, sd=SDdy)^2,
      Phi = rnorm(Npart, object@Synch["Phi"], SDphi)
    )
    attr(df$wFreq, "Synch") <- object@Synch["wFreq"]
    attr(df$wFreq, "SD") <- object@G*object@SD["dy"]^2
    attr(df$Phi, "Synch") <- object@Synch["Phi"]
    attr(df$Phi, "SD") <- object@SD["phi"]
    
    slot(object, "PS") <- df
    
    return(object)
  }
)
setGeneric("fit", def=function(object, sgl) standardGeneric("fit"))
setMethod(
  f="fit", signature = "CSignal",
  definition=function(object){
    n = nrow(object@PS); p0 = object@Synch["Phi"]
    
    f = Sgl ~ n * exp(lam*Time) * sin(w*Time + p0)
    
    nls(f, data=sgl, start=list(lam=-1.4e-3, w=attr(df.p$wFreq, "Synch"))) -> mod1
  }
)