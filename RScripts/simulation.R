source("./RScripts/definitions.R") #this one includes classes.R

library(plyr)
library(ggplot2)
library(mosaic); library(reshape2)

## modeling ##
## changing the signal frequency ####
if(FALSE){
  stu <- CuSampling(Freq=500); mod <- CModel(Phase=pi/32); Ttot = 100
 
  varW0_test(mod, stu, Ttot, wfreqs=c(.05, .3, 1, 3, 5)) -> dat
  
  .stats = ldply(dat, function(e) e$Stats, .id="Freq") #%>% mutate(SE.meanFrq = SE.frq/sqrt(6000*3.6)) 
  
  .stats %>% transmute(Freq, Simulation = SE.frq, Formula = SEAN.frq) %>% 
    reshape2::melt("Freq", variable.name = "Method", value.name = "SE") %>%
    ggplot(aes(Freq, SE, col=Method)) + geom_point() +
    # scale_y_log10() +
    theme_bw() + labs(x=expression(omega), y=expression(sigma[hat(omega)])) + 
    theme(legend.position="top") 
  
  i = 1
  x = dat[[i]]$Sample[seq(1,nrow(dat[[i]]$Sample),length.out=250),]
  ggplot(x, aes(Time, Sgl)) + geom_point() +
    theme_bw() + theme(legend.position="top") + labs(y="signal") +
    geom_line(aes(Time, Sgl), data.frame("Time" = seq(0,Ttot, length.out = 500)) %>% mutate(Sgl = expectation(mod, Time)))
  
  X = ldply(dat, function(e) e$Sample, .id="Freq")
  ggplot(X%>%filter(Freq%in%c("0.05","0.3","3"))) + geom_density(aes(Sgl, col=Freq)) + theme_bw() + labs(x=expression(N[0]*(1+P*exp(lambda*t)*sin(omega*t+phi))))
  
}

## checking the growth of SE with time ####
if(FALSE){
  tau = 100; fs = 500
  
  mod = CModel(decohLam=-1/tau)
  stu = CuSampling(Freq=fs); stm <- CmSampling(Freq=fs, Compaction=.2)
  
  .diagWrap(mod, list("Uni" = stu, "Mod" = stm)) -> dat
  
  dat$Uni$CharPlot + scale_y_log10() + scale_x_log10() + theme_bw()
  dat$Mod$CharPlot + scale_y_log10() + scale_x_log10() + theme_bw()
  
  .stats <- rbind(
    dat$Uni$Stats%>%transmute(Group, SE = SE.frq, SEAN = SEAN.frq, FItot, Smpl="Uni"), 
    dat$Mod$Stats%>%transmute(Group, SE = SE.frq, SEAN = SEAN.frq, FItot, Smpl="Mod")
  ) %>% reshape2::melt(id.vars=c("Smpl","Group","FItot"), variable.name="How", value.name="SE")
  ggplot(.stats, aes(Group, SE, col=Smpl, shape=How)) + geom_point() + theme_bw() + scale_y_log10()
  
}

## number of measurements per observation ####
if(FALSE){
  library(doParallel)
  makeCluster(detectCores()) -> clus; registerDoParallel(clus)
  rtns <- lsf.str(envir=.GlobalEnv, all=TRUE)
  clusterExport(clus, rtns)
  
  tau = 1000; mod = CModel(decohLam=-1/tau)
  stus = llply(c("5" = 5, "50" = 50, "500" = 500, "5000" = 5000), function(f) CuSampling(Freq=f))
  
  ldply(stus, function(stu, .mod, .tau){
    simSample(stu, .mod, .tau) %>% .varWT()
  }, mod, tau, .parallel = TRUE, .paropts = list(.packages="dplyr")) -> .stats; .stats
  # increase in the number of observations has little effect on VarWT so long as the total
  # measurement time remains constant
  # the information Ftr is linearly dependent on the number of observations, confirming the
  # Xtot[zc = s] = n[event/zc] * x0s expression
  
  verr = 10^(0:(nrow(.stats)-1))
  mutate(.stats, VarErr = verr, VarWW = VarErr/Ftr/VarWT) -> .stats; .stats
  # hence, if indeed SD[obs]^2 = SD[meas]^2/n, VarWW should remain the same 
  # and it does
  
  stopCluster(clus)
  
}

## varying the compaction factor while keeping the total time constant ####
if(TRUE){
  library(doParallel)
  makeCluster(detectCores()) -> clus; registerDoParallel(clus)
  rtns <- lsf.str(envir=.GlobalEnv, all=TRUE)
  clusterExport(clus, rtns)
  
  mod = CModel()
  msmpls = llply(c("1.10" = 1.1, ".50" = .5, ".25" = .25, ".10" = .1), 
                 function(cmpt) CmSampling(Freq=500, Compaction=cmpt))
  
  ldply(
    msmpls,
    function(stm, .mod){
      simSample(stm, .mod, 1730) -> smpl
      .varWT(smpl)
    }, mod,
    .parallel = TRUE,
    .paropts = list(.packages = "dplyr"),
    .id = "Compact"
  ) -> .stats; .stats
  
}