source("./RScripts/definitions.R")

library(ggplot2)
library(parallel); library(doParallel); registerDoParallel(detectCores()-1)

## modeling ##
if(FALSE){
  ## 1) trials with different error ####
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
  
  ## 2) checking the growth of omega se with total time ####
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
  
  ## 3) checking the analytical formula ####
  ## the formula is
  ## SEw = SE error / sqrt(sum xj * (sum ti^2 xi/sum xj) - (sum ti xi/sum xj)^2)
  ## this doesn't account for damping
  ## !!!!!!! have to see if damping matters, though !!!!!!
  s = .sample(10); l = nrow(s)
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



