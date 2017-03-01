library(ggplot2)

rm(list=ls(all=TRUE))

source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

lblfnt = 20
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

mod = CModel(); mod@beamLam <- mod@decohLam
smpl = CuSampling()
taud = -1/mod@decohLam
taub = -1/mod@beamLam
tau = taud/(1+taud/taub)
Ttot = 3*tau

.ggplot_Sgl <- function(s){
  ggplot(s[seq(1,nrow(s), length.out=300),]) + 
    geom_line(aes(Time, XSgl), lwd=.2) + 
    geom_point(aes(Time, Sgl), col="red", size=.5) +
    thm +labs(y="Signal")
}

## SIGNAL ####
s = simSample(smpl, mod, Ttot)
.ggplot_Sgl(s)

SNR <- function(model, x, rerr0=3e-2){
  taud = -1/model@decohLam
  taub = -1/model@beamLam
  z = taub/taud
  f = (1-2*z)/2/z
  
  model@Pol/rerr0 * exp(f/taud * x) # here sqrt(dte) from A and from the SNR expression cancel 
                                    # as a result of normalizing A to have rel. error at t=0 rerr0
}

rError <- function(model, x, rerr0=3e-2) rerr0 * exp(-.5*model@beamLam * x) 

## COMPACTION TIME VS W SE ####
varwT <- function(smpl, mod, Ttot){
  sum(smpl$FIDrvt^2) -> denom
  w = smpl$FIDrvt^2/denom
  meanwT = sum(w*smpl$Time)
  
  sum(w*(smpl$Time - meanwT)^2)
}
vwt = varwT(s, mod, Ttot)
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

## SE VS W ####
library(doParallel)
registerDoParallel(detectCores())

llply(c(c(1.57, 3, 6.42)%o%10^(-2:1)), function(w, mod, sampling, duration){
    mod@wFreq <- w
    fit(mod, sampling, duration)
  }, mod, smpl, Ttot, .parallel = TRUE, .paropts = list(.export=lsf.str(), .packages="dplyr")
) -> stats
stopImplicitCluster()

ldply(stats, function(e) e["w", c("Estimate","Std. Error")]) -> stats.w
ggplot(stats.w, aes(Estimate, `Std. Error`)) + geom_point() + 
  thm + labs(x=expression(omega),y=expression(sigma[hat(omega)])) +
  scale_x_log10()


## MODULATION ####
smpl <- CmSampling(sglFreqGuess=rnorm(1,mod@wFreq,1e-4))
simSample(smpl, mod, Ttot) -> s

s <- rbind(s0[,CMPT:="1"], s1[,CMPT:=".33"]%>%dplyr::select(-Node))
ggplot(s[seq(1, nrow(s), length.out = 1300),]) + 
  geom_line(aes(Time, XSgl, col=CMPT), lwd=.1) +
  geom_point(aes(Time, Sgl, col=CMPT), size=.5) +
  thm +
  scale_color_manual(values=c("black","red"))


## 
cmpts <- seq(.1, 1, length.out=5)
names(cmpts) <- as.character(cmpts)
test <- function(cmpt, mod, smpl){
  
  mod@beamLam <- mod@beamLam*as.numeric(cmpt)
  smpl@CMPT <- cmpt
  
  taud = -1/mod@decohLam
  taub = -1/mod@beamLam
  tau = taud/(1+taud/taub)
  
  simSample(smpl, mod, 3*tau) -> s
  s[,CMPT:=as.character(cmpt)]
  #c(Num = nrow(s), Time = s$Time[nrow(s)])
  #fit(mod, smpl, 3*tau)
}

mod = CModel(); mod@beamLam <- mod@decohLam
smpl <- CmSampling(sglFreqGuess=rnorm(1,mod@wFreq,1e-4))

registerDoParallel(detectCores())
ldply(cmpts, function(x, m, s) test(x, m, s), mod, smpl,
      .parallel = TRUE, .paropts = list(.export=lsf.str(), .packages="dplyr")) -> stats
stopImplicitCluster()


stats <- setDT(stats)
stats[,`:=`(Num=length(Time), Last=Time[length(Time)]), by=CMPT]

ggplot(stats[seq(1,nrow(stats), length.out = 2620),]) + geom_line(aes(Time, XSgl, col=CMPT)) +thm

#stats.w <- ldply(stats, function(e) e["w",1:2], .id="CMPT")
#ggplot(stats.w) + geom_pointrange(aes(CMPT, Estimate, ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`))
