library(R6)

RCBunch <- R6Class(
  "RCBunch",
  private = list(G=7e2),
  public = list(
    Synch=c(wFreq=3, Phi=0), SD=numeric(2),
    PS=NULL, Pproj=NULL,
    initialize = function(Npart=1e3, SDdy=1e-3, SDphi=1e-2){
      self$SD <- c(dy = SDdy, phi = SDphi)
      
      self$PS <- data.frame(
        wFreq = self$Synch["wFreq"] + private$G*rnorm(Npart, sd=SDdy)^2,
        Phi = rnorm(Npart, self$Synch["Phi"], SDphi)
      )
      attr(self$PS$wFreq, "Synch") <- self$Synch["wFreq"]
      attr(self$PS$wFreq, "SD") <- private$G*self$SD["dy"]^2
      attr(self$PS$Phi, "Synch") <- self$Synch["Phi"]
      attr(self$PS$Phi, "SD") <- self$SD["phi"]
    },
    Phase = function(at) self$PS$wFreq%o%at + self$PS$Phi,
    project = function(at) self$Pproj <- data.frame(Time=at, Val=colSums( sin(self$Phase(at)) ) )
  ) ## public members
)