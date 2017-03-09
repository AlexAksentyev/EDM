library(R6)
library(data.table)
library(doParallel)
library(dplyr)
library(ggplot2)

RCSignal <- R6Class(
  "RCSignal",
  public = list(
    Bunch=NULL, #reference
    Signal=NULL, specPts=NULL,
    ModelCoef=NULL,
    SDrErr=3e-2,
    initialize=function(bunch, smpl.pts){
      self$Bunch <- bunch
      self$Signal <- self$Bunch$project(smpl.pts)
      private$NSmpl <- 1
      rerr <- rnorm(nrow(self$Signal), sd = self$SDrErr)
      self$Signal[,`:=`(ValNs=Val*(1+rerr), Smpl=private$NSmpl)]
    },
    split=function(){
      
      self$Signal$Time -> smpl.at
      
      registerDoParallel(detectCores())
      
      alply(self$Bunch$EnsPS,1, function(ps, at){
        RCBunch$new(0,NA, NA) -> b
        b$Synch <- unlist(ps); b$EnsPS[, wFreq:=ps[,wFreq]][, Phi:=ps[,Phi]]
        RCSignal$new(b, at)
      }, smpl.at, .parallel = TRUE) -> bl
      
      stopImplicitCluster()
      
      sl
    },
    sample=function(smpl.pts, append=TRUE){
      rerr <- rnorm(length(smpl.pts), sd = self$SDrErr)
      df <- self$Bunch$project(smpl.pts)[,ValNs:=Val*(1+rerr)]
      
      if(append){
        private$NSmpl <- private$NSmpl + 1
        self$Signal <- rbind(
          self$Signal,
          df[,Smpl:=private$NSmpl]
        )
      } else 
        self$Signal <- df[,Smpl:=private$NSmpl]
      
    },
    fit=function(fitpack = NULL){
      if(is.null(self$Signal)) {print("Nothing to fit!"); return(NA)}
      if(is.null(fitpack)) {
        print("Hard-coded")
        n = nrow(self$Bunch$EnsPS); p0 = self$Bunch$Synch["Phi"]; wg = self$Bunch$Synch["wFreq"]
        f = ValNs ~ n * exp(lam*Time) * sin(w*Time + p0)
        guess = list(lam=-1.4e-3, w=wg)
      } else {
        f = fitpack$func; guess = fitpack$guess
      }
      mod3l <- nls(f, data=self$Signal, start=guess)
      self$Signal[,Fit:=fitted(mod3l)]
      coef(summary(mod3l)) -> self$ModelCoef
      return(mod3l)
    },
    findPts=function(what="Node", w.guess=NULL, tol=1e-3){

      k <- switch(what, "Envelope" = TRUE,"Node" = FALSE)
      fn <- function(x) (self$Bunch$project(x)[,Val]^2)
      finder <- function(s, direction, tol){
        dx = pi/self$Bunch$Synch["wFreq"]/2
        i=s[,Time] + c(-dx, +dx)
        optimize(fn, interval=i, maximum = direction, tol=tol)[[1]]->x0
        data.table(self$Bunch$project(x0), Which="Optim")
      }
      
      private$NullSpecPts(what, w.guess) -> pts0
      
      registerDoParallel(detectCores())
      adply(pts0, 1, finder, k, tol,
            .parallel = TRUE, .paropts = list(.packages=c("dplyr","data.table"))) -> pts1
      stopImplicitCluster()
      
      rbind(pts0,pts1) %>% setorder(Time) -> self$specPts
    },
    findNds=function(w.guess=NULL, tol=1e-3){
      private$NullSpecPts(what="Node", w.ref=w.guess) -> pts0
      
      finder <- function(s, tol){
        dx = pi/self$Bunch$Synch["wFreq"]/2
        x0 = s[,Time]
        i =x0 + c(-dx, +dx)
        
        repeat{
          i <- seq(i[1], i[2], length.out=100)
          dx <- i[2]-i[1]
          
          dat = self$Bunch$project(i)
          dat[,`:=`(X0=c(Val[-nrow(dat)], NA), 
                    X1=c(Val[-1], NA))][,Prod:=X0*X1]
          
          x1 = dat[Prod<0,Time]
          
          if(abs(x1-x0) < tol) break
          
          i <- x1 + c(0, dx)
          x0 <- x1
        }
        
        data.table(self$Bunch$project(x1), Which="Optim")
        
      }
      
      adply(pts0,1, finder, tol, .parallel = TRUE, .paropts = list(.packages="data.table")) -> pts1
      
      rbind(pts0,pts1) %>% setorder(Time) -> self$specPts
    },
    Spectrum=function(plot=TRUE){
      require(psd)
      Tstt = self$Signal$Time[1]
      Ttot = self$Signal$Time[nrow(self$Signal)]
      dt = self$Signal$Time[2]-self$Signal$Time[1]
      
      s = ts(self$Signal$ValNs, start=Tstt, end=Ttot, deltat=dt)
      
      pspectrum(s,plot = TRUE) -> sps
      sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)
      
      if(!plot) return(sps)
      
      x = arrange(sps, desc(Pow))[1:20,]
      dw = x[1,"wFreq"]-x[2,"wFreq"]
      
      ggplot(x,aes(wFreq, Pow))+#scale_y_continuous(labels=.fancy_scientific) +
        geom_bar(stat="identity", width=dw*.1) + 
        theme_bw() + labs(x=expression(omega)) +
        geom_vline(xintercept = as.numeric(self$Bunch$Synch["wFreq"]), col="red") -> fps
      
      print(fps)
      return(sps)
    }
  ), ## public
  private = list(
    NSmpl=NULL,
    NullSpecPts=function(what="Node", w.ref=NULL){
      
      w0<- ifelse(!is.null(w.ref), w.ref, self$Bunch$Synch["wFreq"])
      
      p0 = self$Bunch$Synch["Phi"]
      
      d = switch(what, "Envelope" = pi/2, "Node" = 0)
      
      .dum <- function(Time) floor((w0*Time+p0-d)/2/pi)
      
      Nstt = .dum(self$Signal[1, Time])
      Ntot = .dum(self$Signal[nrow(self$Signal), Time])
      
      tnu = (2*pi*Nstt:Ntot-p0+d)/w0; tnu <- tnu[tnu>=0]
      tnd = tnu+pi/w0
      
      data.table(
        N = c(1:length(tnu),1:length(tnd)),
        self$Bunch$project(c(tnu,tnd)), # Time & Val
        Side=rep(c("U","D"),c(length(tnu),length(tnd))),
        Which="Null"
      )
    }
  ) ## private
)

`+.RCSignal` <- function(one, other){
  if(any(one$Bunch$Synch != other$Bunch$Synch)) print("Different synchronous particles!")
  
  at = unique(c(one$Signal$Time, other$Signal$Time))
  
  b0 = RCBunch$new(Npart=0, SDdy=NA, SDphi=NA, WDist="NA")
  b0$Synch <- NULL
  b0$EnsPS <- rbind(one$Bunch$EnsPS, other$Bunch$EnsPS)
  
  RCSignal$new(b0, at)
}
