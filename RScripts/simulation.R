source("./RScripts/definitions.R")

library(plyr)
library(ggplot2)
library(parallel); library(doParallel); n.cores = detectCores()
clus <- makeCluster(n.cores)
registerDoParallel(clus)

## modeling ##
## 1) trials with different error ####
if(FALSE){
  s = .sample(5,15) # sample signal for the two samplings
  # for each type of sampling, and each trial, fit
  s %>% ddply(.(Type, Trl), function(trl) .fit(trl)%>%mutate(SEAN.frq = .compVarF(trl))) -> .stats
  
  .stats%>%group_by(Type) %>% dplyr::summarise(NUM = n(), SE = sd(X.frq), SEAN = SEAN.frq[1], RD = (SE-SEAN)/SEAN)
  
  ## comparing the efficiencies
  daply(.stats, "Type", function(x) .se(x$X.frq)) -> x
  print(x)
  c("SE ratio" = x["uniform"]/x["modulated"], "wg/w0" = w.g/w0) %>% print
  
  ## plot of sample points on signal
  x = filter(s, Trl==1) 
  x %>% mutate(errY = Sgl-XSgl) %>% slice(seq(1,nrow(x), length.out=250)) %>% 
    ggplot(aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
    geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, by=1/fs)) %>% mutate(Sgl=.dcs(Time))) 
}

## 2) checking the growth of omega se with total time ####
if(FALSE){
  s <- .sample(30); nrow(s)/2->l #divide by 2 b/c 2 sampling types
  
  (1:10)*l/10 -> m # split the super-sample into a sequence of cumulatively increasing subsamples;
  # m contains the lengths of the subsamples 
  # for a given length i, pick the first i measurements from each sampling; this is the subsample; combine the samplings
  llply(m, function(i) s%>%ddply("Type", function(x) slice(x,1:i))) -> s
  # fit the subsamples
  ldply(s, function(df) df%>%ddply("Type", .fit)%>%mutate(NUM = nrow(df), SEAN.frq = .compVarF(df))) %>% 
    mutate(SE.MEAN = SE.frq/sqrt(NUM))->.stats
  
  ggplot(.stats, aes(NUM, SE.MEAN, col=Type)) + geom_point() + 
    scale_y_log10() +
    theme_bw() + labs(x="Sample size",y=expression(sigma[bar(omega)])) + theme(legend.position="top")
  
  ggplot(mutate(.stats, Dev = (SE.frq-SEAN.frq)/SEAN.frq), aes(NUM, abs(Dev), col=Type)) + geom_point() + 
    theme_bw() + labs(x="Sample size", y="(Simulation - Analytic)/Analytic SE") + theme(legend.position="top") +
    scale_y_log10()
  
  ggplot(.stats, aes(SEAN.frq, SE.frq, col=Type)) + geom_point() + 
    theme_bw() + labs(x="Analytic SE", y="Simulation SE") + theme(legend.position="top") + 
    geom_abline(slope=1,intercept=0) + 
    scale_y_log10() + scale_x_log10()
}

## 3) checking the analytical formula ####
if(FALSE){
  ## the formula is
  ## SEw = SE error / sqrt(sum xj * (sum ti^2 xi/sum xj) - (sum ti xi/sum xj)^2)
  ## this doesn't account for damping
  ## !!!!!!! have to see if damping matters, though !!!!!!
  s = .sample(480); #l = nrow(s)
  .stats <- ddply(s,"Type", .fit, .parallel = TRUE) %>% mutate(SEAN.frq = daply(s, "Type", .compVarF))
  paste("Compaction factor: ", comptn) %>% print()
  .stats%>% print()
  paste("SE ratio uni/mod", round(.stats[1,"SE.frq"]/.stats[2,"SE.frq"],2)) %>% print()
  paste("SEAN ratio uni/mod", round(.stats[1,"SEAN.frq"]/.stats[2,"SEAN.frq"],2)) %>% print()
  paste("SE/SEAN ratio, uni", round(.stats[1,"SE.frq"]/.stats[1,"SEAN.frq"],2))%>%print()
  paste("SE/SEAN ratio, mod", round(.stats[2,"SE.frq"]/.stats[2,"SEAN.frq"],2))%>%print()
  
  x = mutate(s, errY = Sgl-XSgl) %>% slice(seq(1,l, length.out=250))
  ggplot(x, aes(Time, Sgl)) + geom_pointrange(aes(ymin=Sgl-errS, ymax=Sgl+errS,col=Type), size=.3) + theme_bw() + 
    geom_line(aes(Time, Sgl), data.frame(Time = seq(0, Ttot, length.out=250)) %>% mutate(Sgl=.dcs(Time))) +
    theme(legend.position="top") ->p
  p
}

## 4) testing different compaction factors ####
if(FALSE){
  c(seq(.5,.1,-.1), seq(.05,.01,-.01)) -> alphas;# names(alphas) <- alphas
  alphas = c(1:10)%o%10^(-2:-1) %>%c
  ldply(alphas, function(al) { 
    dat=.msampleF(Nprd=45, comptn = al, len=15600);
    data.frame("SE" = .fit(dat)["SE.frq"], "SEAN.frq" = .compVarF(dat))
  },.parallel=TRUE) %>% cbind("Comptn"=alphas) %>% melt(id.vars="Comptn", variable.name="Which") -> x
  ggplot(x, aes(Comptn,value, col=Which)) + geom_point() + 
    scale_y_log10() + scale_x_log10() +
    theme_bw() + labs(x="compaction factor", y=expression(sigma[hat(omega)])) + 
    scale_color_discrete(breaks=c("SE.frq","SEAN.frq"), labels=c("Simulation","Formula")) -> p
  
  p
}

## 5) changing the signal frequency ####
if(TRUE){
  stu <- usampling(freq=50); Ttot = 100
  varW0_test <- function(.mod, .Time){
    # cat(paste("model frequency", .mod$wfreq, "\n"))
    
    sample.usampling(how = stu, model = .mod, how.long = .Time) -> .spl
    # cat(paste("sample size", nrow(.spl), "\n"))
    
    .fit(.spl, .mod) -> .stats
    # .stats = NULL
    
    list("Stats" = .stats, "Sample" = .spl)
  }
  
  mod <- model(phs=pi/2); 
  w0s = mod$wfreq*c(.01, .1,.5, 1, 5, 10); names(w0s) <- w0s
  llply(w0s, function(w) setWFreq(mod, w)) -> mods
  rtns <- c("stu",".extract.stats", 
            "sample","sample.usampling",
            "derivative", "derivative.model",
            "expectation", "expectation.model",
            ".fit", "varW0_test"
           )
  # clusterExport(clus, rtns)
  llply(
    mods, varW0_test, Ttot, 
    .inform = FALSE,
    .parallel=TRUE, 
    .paropts = list(.packages=c("dplyr"), .export = rtns)
  ) -> dat
  
  .stats = ldply(dat, function(e) e$Stats, .id="Freq"); .stats
  
  ggplot(.stats, aes(Freq, SE.frq)) + geom_point() + 
    scale_y_log10() +
    theme_bw() + labs(x=expression(omega), y=expression(sigma[hat(omega)])) + 
    theme(legend.position="top") 
  
  i = "15"
  x = dat[[i]]$Sample[seq(1,nrow(dat[[i]]$Sample),length.out=250),]
  ggplot(x, aes(Time, Sgl)) + geom_point() +
    theme_bw() + theme(legend.position="top") + labs(y="signal") +
    geom_line(aes(Time, Sgl), data.frame("Time" = seq(0,Ttot, length.out = 500)) %>% mutate(Sgl = expectation(mod, Time)))
  
  X = ldply(dat, function(e) e$Sample, .id="Freq")
  ggplot(X%>%filter(Freq%in%c("15","30"))) + geom_density(aes(Sgl, col=Freq)) + theme_bw() + labs(x=expression(N[0]*(1+P*exp(lambda*t)*sin(omega*t+phi))))
  
}


## 
