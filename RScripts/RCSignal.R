library(R6)
library(data.table)
library(doParallel)
library(dplyr)

RCSignal <- R6Class(
  "RCSignal",
  public = list(
    Bunch=NULL, #reference
    Signal=NULL, specPts=NULL,
    ModelCoef=NULL,
    initialize=function(bunch, smpl.pts){
      self$Bunch <- bunch
      self$Signal <- self$Bunch$project(smpl.pts)
      private$NSmpl <- 1
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
      private$NSmpl <- private$NSmpl + 1
      self$Signal <- rbind(
        self$Signal[,Smpl:=private$NSmpl-1],
        self$Bunch$project(smpl.pts)[,Smpl:=private$NSmpl]
      )
    },
    fit=function(fitpack = NULL){
      if(is.null(self$Signal)) {print("Nothing to fit!"); return(NA)}
      if(is.null(fitpack)) {
        print("Hard-coded")
        n = nrow(self$Bunch$EnsPS); p0 = self$Bunch$Synch["Phi"]; wg = self$Bunch$Synch["wFreq"]
        f = Val ~ n * exp(lam*Time) * sin(w*Time + p0)
        guess = list(lam=-1.4e-3, w=wg)
      } else {
        f = fitpack$func; guess = fitpack$guess
      }
      coef(summary(nls(f, data=self$Signal, start=guess))) -> self$ModelCoef
      return(self$ModelCoef)
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
      
      rbind(pts0,pts1) %>% arrange(Time) -> self$specPts
    },
    Spectrum=function(plot=TRUE, method="ar"){
      Tstt = self$Signal$Time[1]
      Ttot = self$Signal$Time[nrow(self$Signal)]
      dt = self$Signal$Time[2]-self$Signal$Time[1]
      
      s = ts(self$Signal$Val, start=Tstt, end=Ttot, deltat=dt)
      
      method <- eval(parse(text=paste0("spec.", method)))
      method(s,plot = FALSE) -> sps
      sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)
      
      if(!plot) return(sps)
      
      x = arrange(sps, desc(Pow))[1:20,]
      dw = x[1,"wFreq"]-x[2,"wFreq"]
      
      ggplot(x,aes(wFreq, Pow))+#scale_y_continuous(labels=.fancy_scientific) +
        geom_bar(stat="identity", width=dw*.1) + 
        theme_bw() + labs(x=expression(omega)) +
        geom_vline(xintercept = as.numeric(self$Bunch$Synch["wFreq"]), col="red") -> fps
      
      fps
      return(sps)
    }
  ), ## public
  private = list(
    NSmpl=NULL,
    NullSpecPts=function(what="Node", w.ref=NULL){
      
      w0<- ifelse(!is.null(w.ref), w.ref, self$Bunch$Synch["wFreq"])
      
      p0 = self$Bunch$Synch["Phi"]
      
      .dum <- function(Time) floor((w0*Time+p0-pi/2)/2/pi)
      
      Nstt = .dum(self$Signal[1, Time])
      Ntot = .dum(self$Signal[nrow(self$Signal), Time])
      
      d = switch(what, "Envelope" = pi/2, "Node" = 0)
      
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
