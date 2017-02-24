library(R6)
library(data.table)

RCBunch <- R6Class(
  "RCBunch",
  private = list(
    G=7e2, SD=numeric(2),
    polProj=function(at){
      n = nrow(self$EnsPS)
      # do computation one measurement at a time to minimize memory use
      registerDoParallel(detectCores())
      aaply(at,1, function(x, npart){ sum(sin(self$EnsPS[,"wFreq"]*x+self$EnsPS[,"Phi"]))
        ## checked near equivalence => no need for it now
        # r = runif(npart, 0, 1e-1)
        # c(
        #   Val = sum(sin(self$EnsPS[,"wFreq"]*x+self$EnsPS[,"Phi"])),
        #   Valu = sum(sin(self$EnsPS[,"wFreq"]*x+self$EnsPS[,"Phi"]+ self$EnsPS[,"wFreq"]*r))
        # )
      }, n, .parallel = TRUE) ->res
      stopImplicitCluster()
      
      res
    }
  ), ## private members
  public = list(
    Synch=c(wFreq=3, Phi=0),
    EnsPS=NULL, # keeps ensemble phase space
    initialize = function(Npart=1e3, SDdy=1e-3, SDphi=1e-2, WDist="phys", ...){
      require(bigmemory)
      # if(Npart>1e4) {print("Will have computational problems; truncating."); Npart <- 1e4}
      private$SD <- c(dy = SDdy, phi = SDphi)
      
      # to get the same scale of chisq and norm distributions
      q0 = qchisq(pnorm(1)-pnorm(-1),1) #the quantile in chisq corresponding to 1 sd norm
      # dw/(G*sddy^2) = q0 => dw = q0*G*sddy^2
      # dw(1sd) = sdw, and since norm is symmetric, sdw ~ q0/2
      
      supplied = list(...); sname = names(supplied)
      sdw = ifelse("SDwFreq" %in% sname, supplied$SDwFreq, .5*q0*private$G*SDdy^2)
      shpw = ifelse("ShpwFreq" %in% sname, supplied$ShpFreq, 1.5)
      skeww = ifelse("SkewwFreq" %in% sname, supplied$SkewwFreq, 1)
      
      self$EnsPS <- big.matrix(nrow=Npart, ncol=2, dimnames = list(NULL, c("wFreq","Phi")))
      
      self$EnsPS[,"wFreq"] <- switch(WDist,
        "phys" = self$Synch["wFreq"] + private$G*rnorm(Npart, sd=SDdy)^2,
        "norm" = rnorm(Npart, self$Synch["wFreq"], sdw),
        "skew" = self$Synch["wFreq"] + skeww*(
          rweibull(Npart, shpw, sdw) - 
            ifelse(shpw>1, sdw*((shpw-1)/shpw)^(1/shpw), 0)
        ),
        numeric(0)
      )
      
      self$EnsPS[,"Phi"] <- rnorm(Npart, self$Synch["Phi"], SDphi)
      
      # self$EnsPS <- data.table(wFreq=w, Phi=rnorm(Npart, self$Synch["Phi"], SDphi))
      # 
      # attr(self$EnsPS$wFreq, "Synch") <- self$Synch["wFreq"]
      # attr(self$EnsPS$wFreq, "SD") <- private$G*private$SD["dy"]^2
      # attr(self$EnsPS$Phi, "Synch") <- self$Synch["Phi"]
      # attr(self$EnsPS$Phi, "SD") <- private$SD["phi"]
    },
    Phase = function(at) self$EnsPS[,"wFreq"]%o%at + self$EnsPS[,"Phi"],
    project = function(at) data.table(Time=at, Val=private$polProj(at) )
  ) ## public members
)

