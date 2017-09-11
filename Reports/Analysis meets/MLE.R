library(plyr); library(dplyr); library(ggplot2); library(reshape2)
library(data.table)

w0 = 3; se = 7e-2
w = w0*c(1.01,.95,1,1.05); names(w) <- paste0("S",0:(length(w)-1))
Time = seq(0,6,.01)

llply(w,function(h) sin(h*Time)) -> s
do.call(cbind,s) %>% cbind(Time) %>% data.table() -> s
s[, S0 := S0 + rnorm(nrow(s),0,se)]
melt(s, id.vars="Time", variable.name = "model") -> s

ggplot(s[model!="S0"], aes(Time, value, col=model)) + geom_line() + 
  geom_pointrange(aes(Time, value, ymin=value-se,ymax=value+se), col="black", 
                  data = s[model=="S0"][seq(1,length(Time),length.out = 10),]) +
  theme_bw() + theme(legend.position="top")
