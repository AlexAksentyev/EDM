library(R6)
library(bigmemory)

# should make bunch$Phase() return a bigmemory matrix

RCBunch <- R6Class(
  "RCBunch",
  private = list(
    G=7e2, SD=numeric(2),
    # polProj=function(at) private$big_aggregate(self$Phase(at), function(X) colSums(sin(X)), .combine = 'c'),
    polProj=function(at) colSums(sin(self$Phase(at)[,])),
    # utilities
    cutBySize=function(m, block.size, nb = ceiling(m / block.size)) {
      int <- m / nb
      
      upper <- round(1:nb * int)
      lower <- c(1, upper[-nb] + 1)
      size <- c(upper[1], diff(upper))
      
      cbind(lower, upper, size)
    },
    seq2=function(lims) seq(lims["lower"], lims["upper"]),
    big_aggregate=function(X, FUN, .combine, block.size = 1e3) {
      require(foreach)
      intervals <- private$cutBySize(ncol(X), block.size)
      
      foreach(k = 1:nrow(intervals), .combine = .combine) %do% {
        FUN(X[, private$seq2(intervals[k, ])])
      }
    }
  ), ## private members
  public = list(
    Synch=c(wFreq=3, Phi=0),
    EnsPS=NULL, # keeps ensemble phase space
    initialize = function(Npart=1e3, SDdy=1e-3, SDphi=1e-2, WDist="phys", ...){
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
    Phase = function(at) as.big.matrix(self$EnsPS$wFreq%o%at + self$EnsPS$Phi),
    project = function(at) data.table(Time=at, Val=private$polProj(at) )
  ) ## public members
)