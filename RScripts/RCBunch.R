library(R6)
library(ggplot2)
library(dplyr)

RCBunch <- R6Class(
  "RCBunch",
  private = list(
    G=7e2,
    Func=function(at) colSums( sin(self$Phase(at)) ),
    NullSpecPts=function(what="Node", w.ref=NULL){
      
      if(!is.null(w.ref)) w0 <- w.ref
      else w0 <- self$Synch["wFreq"]
      
      p0 = self$Synch["Phi"]
      
      .dum <- function(Time) floor((w0*Time+p0-pi/2)/2/pi)
      
      Nstt = .dum(self$Pproj[1, "Time"])
      Ntot = .dum(self$Pproj[nrow(self$Pproj), "Time"])
      
      d = switch(what, "Envelope" = pi/2, "Node" = 0)
      
      tnu = (2*pi*Nstt:Ntot-p0+d)/w0; tnu <- tnu[tnu>=0]
      tnd = tnu+pi/w0
      
      data.frame(
        N = c(1:length(tnu),1:length(tnd)),
        Time=c(tnu, tnd),
        Val=private$Func(c(tnu,tnd)), 
        Side=rep(c("U","D"),c(length(tnu),length(tnd))),
        Which="Null"
      )
    }
  ), ## private
  public = list(
    Synch=c(wFreq=3, Phi=0), SD=numeric(2),
    PS=NULL, Pproj=NULL, specPts=NULL,
    Model=NULL,
    initialize = function(Npart=1e3, SDdy=1e-3, SDphi=1e-2, WDist="phys", ...){
      self$SD <- c(dy = SDdy, phi = SDphi)
      
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
        )
      )
      
      self$PS <- data.frame(
        wFreq = w,
        Phi = rnorm(Npart, self$Synch["Phi"], SDphi)
      )
      attr(self$PS$wFreq, "Synch") <- self$Synch["wFreq"]
      attr(self$PS$wFreq, "SD") <- private$G*self$SD["dy"]^2
      attr(self$PS$Phi, "Synch") <- self$Synch["Phi"]
      attr(self$PS$Phi, "SD") <- self$SD["phi"]
    },
    Phase = function(at) self$PS$wFreq%o%at + self$PS$Phi,
    project = function(at) self$Pproj <- data.frame(Time=at, Val=private$Func(at) ),
    fit = function(fitpack = NULL){
      if(is.null(self$Pproj)) {print("Nothing to fit!"); return(NA)}
      if(is.null(fitpack)) {
        print("Hard-coded")
        n = nrow(self$PS); p0 = self$Synch["Phi"]; wg = self$Synch["wFreq"]
        f = Val ~ n * exp(lam*Time) * sin(w*Time + p0)
        guess = list(lam=-1.4e-3, w=wg)
      } else {
        f = fitpack$func; guess = fitpack$guess
      }
      nls(f, data=self$Pproj, start=guess) -> self$Model
      return(summary(self$Model))
    },
    findPts=function(what="Node", w.guess=NULL, tol=1e-3){
      require(doParallel); n.cores = detectCores()
      clus <- makeCluster(n.cores)
      registerDoParallel(clus)
      
      k <- switch(what, "Envelope" = TRUE,"Node" = FALSE)
      fn <- function(x) (private$Func(x))^2
      finder <- function(s, direction, tol){
        dx = pi/self$Synch["wFreq"]/2
        i=as.numeric(c(s["Time"]-dx, s["Time"]+dx))
        optimize(fn, interval=i, maximum = direction, tol=tol)[[1]]->x0
        data.frame(Time=x0, Val=private$Func(x0), Which="Optim")
      }
      
      private$NullSpecPts(what, w.guess) -> pts0
      adply(pts0, 1, finder, k, tol,
            .parallel = FALSE, .paropts = list(.packages="dplyr")) -> pts1
      
      stopCluster(clus)

      rbind(pts0,pts1) %>% arrange(Time) -> self$specPts
    },
    Spectrum=function(plot=TRUE, method="ar"){
      Tstt = self$Pproj$Time[1]
      Ttot = self$Pproj$Time[nrow(self$Pproj)]
      dt = self$Pproj$Time[2]-self$Pproj$Time[1]
      
      s = ts(self$Pproj$Val, start=Tstt, end=Ttot, deltat=dt)
      
      method <- eval(parse(text=paste0("spec.", method)))
      method(s,plot = FALSE) -> sps
      sps <- data.frame(Freq=sps$freq, Pow=sps$spec) %>% mutate(wFreq=2*pi*Freq)
      
      if(!plot) return(sps)
      
      x = arrange(sps, desc(Pow))[1:20,]
      dw = x[1,"wFreq"]-x[2,"wFreq"]
      
      sdw = self$SD["wFreq"]
      
      ggplot(x,aes(wFreq, Pow))+#scale_y_continuous(labels=.fancy_scientific) +
        geom_bar(stat="identity", width=dw*.1) + 
        theme_bw() + labs(x=expression(omega)) +
        geom_vline(xintercept = self$Synch["wFreq"], col="red") %>% print
      return(sps)
    }
  ) ## public members
)