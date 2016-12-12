source("./RScripts/definitions.R") #this one includes classes.R

library(plyr)
library(ggplot2)

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
if(TRUE){
  library(mosaic) #for derivedFactor
  tau = 100; fs = 500
  
  mod = CModel(decohLam=-1/tau)
  stu = CuSampling(Freq=fs); stm <- CmSampling(Freq=fs, Compaction=.2)
  
  ## this here will go into a single test function akin to varW0_test ##
  ## so, give it a model and maybe a list of samplings, and then compute the samples,
  ## split them, and do the analysis
  
  ## .diagnosis code in here as well
  # 
  # simSample(stu, mod, Ttot, rerror=err) %>% mutate(Group = derivedFactor(
  #   "A" = Time <= sptt[1],
  #   "B" = Time <= sptt[2],
  #   "C" = Time <= sptt[3],
  #   "D" = Time <= sptt[4],
  #   "E" = Time <= sptt[5],
  #   .method = "first"
  # )) -> usmpl
  # 
  # simSample(stm, mod, Ttot, rerror=err) %>% mutate(Group = derivedFactor(
  #   "A" = Time <= sptt[1],
  #   "B" = Time <= sptt[2],
  #   "C" = Time <= sptt[3],
  #   "D" = Time <= sptt[4],
  #   "E" = Time <= sptt[5],
  #   .method = "first"
  # )) -> msmpl
  # 
  # udat <- .diagnosis(usmpl)
  # mdat <- .diagnosis(msmpl)
  
  ## and from here on it's the results, so it'll remain in simulation.R ##

  .diagWrap(mod, list("Uni" = stu, "Mod" = stm)) -> dat
  
  dat$Uni$CharPlot + scale_y_log10() + scale_x_log10() + theme_bw()
  dat$Mod$CharPlot + scale_y_log10() + scale_x_log10() + theme_bw()
  
  .stats <- rbind(
    dat$Uni$Stats%>%transmute(Group, SE = SE.frq, SEAN = SEAN.frq, FItot, Smpl="Uni"), 
    dat$Mod$Stats%>%transmute(Group, SE = SE.frq, SEAN = SEAN.frq, FItot, Smpl="Mod")
  ) %>% reshape2::melt(id.vars=c("Smpl","Group","FItot"), variable.name="How", value.name="SE")
  ggplot(.stats, aes(Group, SE, col=Smpl, shape=How)) + geom_point() + theme_bw()
  
}

## 
