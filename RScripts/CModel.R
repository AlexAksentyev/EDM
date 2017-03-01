library(methods)
library(dplyr); library(plyr)

CModel = setClass(
  "CModel",
  slots = c(wFreq="numeric", Phase="numeric", Num0="numeric", beamLam = "numeric", Pol="numeric", decohLam="numeric"),
  prototype = list(wFreq=3, Phase=0, Num0=6730, beamLam=-1/2000, Pol=.4, decohLam=log(.25)/1000)
)

setGeneric("setValue",def=function(object, value) standardGeneric("setValue"))
setGeneric("expectation", def=function(object, at) standardGeneric("expectation"))
setGeneric("fiDer", def=function(object, at) standardGeneric("fiDer"))
setGeneric("timeDer", def=function(object, at) standardGeneric("timeDer"))
setGeneric("nodes", def=function(object, number) standardGeneric("nodes"))
setGeneric("fit", def=function(object, sampling, duration) standardGeneric("fit"))

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
  f="nodes", signature = "CModel",
  definition = function(object, number){
    w = object@wFreq; phi = object@Phase
    
    (pi*0:number-phi)/w
  }
)
setMethod(
  f="fit", signature = "CModel",
  definition = function(object, sampling, duration){
    
    N0 <- object@Num0; P <- object@Pol; w0 <- object@wFreq; phi <- object@Phase
    lam.decoh <- object@decohLam; lam.beam <- object@beamLam
    
    f = Sgl ~ N0 * exp(lb*Time) * (1 + P*exp(ld*Time)*sin(w*Time + phi))
    guess = llply(c("w" = w0, "lb"= lam.beam, "ld" = lam.decoh), function(x) rnorm(1, x, abs(x)*1e-4))
    
    s = simSample(sampling, object, duration)
    nls(f, s, guess) %>% summary %>% coef
  }
)