library(R6)
library(data.table)

RCBunch <- R6Class(
  "RCBunch",
  private = list(
    G=7e2, SD=numeric(2),
    polProj=function(at) colSums(sin(self$Phase(at)))
  ), ## private members
  public = list(
    Synch=c(wFreq=3, Phi=0),
    EnsPS=NULL, # keeps ensemble phase space
    initialize = function(Npart=1e3, SDdy=1e-3, SDphi=1e-2, WDist="phys", ...){
      if(Npart>5e4) {print("Will have computational problems; truncating."); Npart<-5e4}
      private$SD <- c(dy = SDdy, phi = SDphi)
      
      supplied = list(...); sname = names(supplied)
      sdw = ifelse("SDwFreq" %in% sname, supplied$SDwFreq, private$G*SDdy^2)
      shpw = ifelse("ShpwFreq" %in% sname, supplied$ShpFreq, 1.5)
      skeww = ifelse("SkewwFreq" %in% sname, supplied$SkewwFreq, 1)
      
      w = switch(WDist,
        "phys" = self$Synch["wFreq"] + private$G*rnorm(Npart, sd=SDdy)^2,
        "norm" = rnorm(Npart, self$Synch["wFreq"], sdw),
        "skew" = self$Synch["wFreq"] + skeww*(
          rweibull(Npart, shpw, sdw) - 
            ifelse(shpw>1, sdw*((shpw-1)/shpw)^(1/shpw), 0)
        ),
        numeric(0)
      )
      
      self$EnsPS <- data.table(
        wFreq = w,
        Phi = rnorm(Npart, self$Synch["Phi"], SDphi)
      )
      attr(self$EnsPS$wFreq, "Synch") <- self$Synch["wFreq"]
      attr(self$EnsPS$wFreq, "SD") <- private$G*private$SD["dy"]^2
      attr(self$EnsPS$Phi, "Synch") <- self$Synch["Phi"]
      attr(self$EnsPS$Phi, "SD") <- private$SD["phi"]
    },
    Phase = function(at) self$EnsPS$wFreq%o%at + self$EnsPS$Phi%o%rep(1,length(at)),
    project = function(at) data.table(Time=at, Val=private$polProj(at) )
  ) ## public members
)