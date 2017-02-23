library(ggplot2)

source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

lblfnt = 20
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

mod = CModel()
smpl = CuSampling()
Ttot = -2.3/mod@decohLam

varwT <- function(smpl, mod, Ttot){
  simSample(smpl, mod, Ttot) -> s
  sum(s$FIDrvt^2) -> denom
  w = s$FIDrvt^2/denom
  meanwT = sum(w*s$Time)
  
  sum(w*(s$Time - meanwT)^2)
}
vwt = varwT(smpl, mod, Ttot)
g <- function(x) {
  lam = mod@decohLam
  lamp = mod@decohLam*pi/mod@wFreq
  
  ifelse(lam != 0, (exp(lam*x)-1)/(exp(lamp)-1), x*mod@wFreq/pi)
}
Gtot = g(Ttot)

SDe = smpl@rerror*mod@Num0

dte = 1e-5
SEw = 1e-4
Xtot = (SDe/SEw)^2/vwt

# f <- function(dt){
#   w = mod@wFreq
#   dt + sin(w*dt)/w - 2*Xtot*dte/Gtot
# }
# 
# uniroot(f,c(0,1))$root -> x0

# x = seq(0, 1, 1e-3)
# plot(f(x)~x, type="l"); abline(h=0, v=x0, col="red"); abline(v=pi/mod@wFreq/2, lty=2, col="gray")

as.numeric(c(1,2,5)%o%10^(-3:-6))->sew
ldply(sew,
  function(sew){
    xtot = (SDe/sew)^2/vwt
    .f <- function(dt){
      w = mod@wFreq
      dt + sin(w*dt)/w - 2*xtot*dte/Gtot
    }
    
    x0 <- tryCatch(uniroot(.f, c(0,1))$root, error = function(e) x0 <- NA)
    
    c("SEw" = sew,"Xtot" = xtot, "Dt" = x0, "Dt_t" = x0/pi*mod@wFreq*100)
  }, .id=NULL
) -> tdf

tdf.p <- filter(tdf, !is.na(Dt)) %>%arrange(SEw)
tdf.p <- transmute(tdf.p, SEw = factor(SEw), CMPT = Dt, Lev = 1-seq(1,0, length.out=nrow(tdf.p)))
  
at = seq(-.5, 2, length.out = 100)
expectation(mod, at) -> sgl
sdf = data.frame(Time=at, Val=sgl)
prd <- pi/mod@wFreq*c(.5, 1)
ggplot(sdf, aes(Time, Val)) + geom_line(lwd=.5) + thm +
  scale_color_discrete(name=expression(sigma[hat(omega)])) +
  labs(y="Signal") +
  geom_hline(yintercept = mod@Num0, lty=3) + 
  geom_vline(xintercept=c(0, prd), lty=3) +
  geom_segment(aes(x=0, xend=CMPT, y=mod@Num0-1500*Lev, yend=mod@Num0-1500*Lev, col=SEw), arrow = arrow(ends="last",length = unit(.05,"inches"), angle=90), data=tdf.p) +
  geom_segment(x=0, xend=prd[1], y=mod@Num0*.73, yend=mod@Num0*.73, arrow = arrow(ends="last", angle=90, length=unit(.1,"inches"))) + 
  annotate("text", x=prd[1]/2, y=mod@Num0*.75, label="Uniform sampling")
