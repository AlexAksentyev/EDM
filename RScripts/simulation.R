source("./RScripts/classes.R")
source("./RScripts/definitions.R")

library(plyr)
library(ggplot2)

## modeling ##
## changing the signal frequency ####
if(FALSE){
  stu <- CuSampling(Freq=5000); mod <- CModel(Phase=0); Ttot = 100
 
  varW0_test(mod, stu, Ttot, wfreqs=c(.05, .3, 1, 3, 5)) -> dat
  
  .stats = ldply(dat, function(e) e$Stats, .id="Freq") #%>% mutate(SE.meanFrq = SE.frq/sqrt(6000*3.6)) 
  
  .stats %>% transmute(Freq, Simulation = SE.frq, Formula = SEAN.frq) %>% 
    reshape2::melt("Freq", variable.name = "Method", value.name = "SE") %>%
    ggplot(aes(Freq, SE, col=Method)) + geom_point() +
    # scale_y_log10() +
    theme_bw() + labs(x=expression(omega), y=expression(sigma[hat(omega)])) + 
    theme(legend.position="top") 
  
  i = "1"
  x = dat[[i]]$Sample[seq(1,nrow(dat[[i]]$Sample),length.out=250),]
  ggplot(x, aes(Time, Sgl)) + geom_point() +
    theme_bw() + theme(legend.position="top") + labs(y="signal") +
    geom_line(aes(Time, Sgl), data.frame("Time" = seq(0,Ttot, length.out = 500)) %>% mutate(Sgl = expectation(mod, Time)))
  
  X = ldply(dat, function(e) e$Sample, .id="Freq")
  ggplot(X%>%filter(Freq%in%c("0.01","0.3","3"))) + geom_density(aes(Sgl, col=Freq)) + theme_bw() + labs(x=expression(N[0]*(1+P*exp(lambda*t)*sin(omega*t+phi))))
  
}

## checking the growth of SE with time ####
if(TRUE){
  library(mosaic)
  Ttot = 500; fs = 500; err = 1e-1
  
  mod = CModel(decohLT=log(.25)/Ttot)
  stu = CuSampling(Freq=fs); stm <- CmSampling(Freq=fs, Compaction=.2)
  
  .diagnosis <- function(smpl){
    smpl <- mutate(smpl, Char = abs(Drvt)/(mod@Num0*mod@Pol))
    smplslc <- smpl[seq(1,nrow(smpl), length.out=500),]
    psigscat <- ggplot(smplslc, aes(Time, Sgl, col=Group)) + geom_point() + 
      geom_line(
        aes(Time, XSgl), linetype=3, 
        data=data.frame(Group=NA, Time = seq(0, smpl[nrow(smpl), "Time"], length.out=500))%>%
          mutate(XSgl = expectation(mod, Time))
      ) + theme_bw()
    
    ddply(smpl, "Group", function(g) {.fit(g, mod)%>%mutate(SEAN.frq = .compVarF(smpl), Fish = mean(g$Drvt^2), Smpl="U")})  -> .stats
    .stats %>% transmute(Group, SE = SE.frq) -> .stats
    pse <- ggplot(.stats, aes(Group, SE)) + geom_point() + theme_bw()
    
    pchg <- ggplot(smpl) + geom_density(aes(Char, col=Group))
    
    list("SglPlot" = psigscat, "SEPlot" = pse, "CharPlot" = pchg, "Sample" = smpl, "Stats" = .stats)
  }
  
  sptt <- seq(.25,1,.25)*Ttot
  
  simSample(stu, mod, Ttot, rerror=err) %>% mutate(Group = derivedFactor(
    "A" = Time <= sptt[1],
    "B" = Time <= sptt[2],
    "C" = Time <= sptt[3],
    "D" = Time <= sptt[4],
    .method = "first"
  )) -> usmpl
  
  
  simSample(stm, mod, Ttot, rerror=err) %>% mutate(Group = derivedFactor(
    "A" = Time <= sptt[1],
    "B" = Time <= sptt[2],
    "C" = Time <= sptt[3],
    "D" = Time <= sptt[4],
    .method = "first"
  )) -> msmpl
  
  udat <- .diagnosis(usmpl)
  mdat <- .diagnosis(msmpl)
  
  .stats <- rbind(udat$Stats%>%mutate(Smpl="Uni"), mdat$Stats%>%mutate(Smpl="Mod"))
  ggplot(.stats, aes(Group, SE, col=Smpl)) + geom_point() + theme_bw()
}

## 
