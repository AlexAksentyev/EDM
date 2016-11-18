library(plyr); library(dplyr); library(reshape2)
library(ggplot2)
library(Hmisc)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

Prd = 1000
st = seq(3,Prd, by=.01)
f = c(1,2,5)%o%10^(c(-5))%>%c; names(f) <- as.character(f)

wt = ldply(st, function(t) (cos(f*t + pi/4))^2) %>% 
  melt(variable.name="Freq", value.name="X") %>%
  mutate(Cat = cut(X, breaks=2, labels=c("0","1"))) %>% 
  cbind(data.frame("STime" = rep(st, times=length(f)))) %>%
  ddply("Freq", function(s) mutate(s, Wt = X/sum(X)))

ggplot(wt, aes(X, col = Freq)) + geom_density() + 
  scale_x_continuous(name=expression(cos^2~(omega%.%t), ")")) +
  theme_minimal() +
  ggtitle("Point Fisher information distributions for different signal frequencies")


nu = 1e-3
ddply(wt, "Freq", function(s){
  wt.m.t = with(s, wtd.mean(STime, Wt))
  wt.v.t = with(s, sum(Wt*(STime - wt.m.t)^2))
  ftr = sum(s$X)
  c("Factor" = ftr, "Wt.Mean.t" = wt.m.t, "Wt.Var.t" = wt.v.t,"varW" = nu/(ftr*wt.v.t))
})-> varW.w

ggplot(varW.w, aes(Freq, varW/min(varW))) + geom_point() + theme_minimal() +
  scale_y_continuous(name=expression(sigma^2~group("(",list(hat(omega), omega),")"))) +
  scale_x_discrete(name=expression(omega), labels = fancy_scientific) + 
  ggtitle("Variance of the omega estimate conditional on the value") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

