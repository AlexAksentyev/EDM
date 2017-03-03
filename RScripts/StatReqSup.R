library(ggplot2)
library(doParallel)
library(reshape2)

rm(list=ls(all=TRUE))

source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

lblfnt = 20
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

varwT <- function(smpl, mod){
  sum(smpl$FIDrvt^2) -> denom
  w = smpl$FIDrvt^2/denom
  meanwT = sum(w*smpl$Time)
  
  sum(w*(smpl$Time - meanwT)^2)
}
g <- function(x) {
  lam = mod@decohLam
  lamp = mod@decohLam*pi/mod@wFreq
  
  ifelse(lam != 0, (exp(lam*x)-1)/(exp(lamp)-1), x*mod@wFreq/pi)
}
SNR <- function(model, x, rerr0=3e-2){
  taud = -1/model@decohLam
  taub = -1/model@beamLam
  z = taub/taud
  f = (1-2*z)/2/z
  
  model@Pol/rerr0 * exp(f/taud * x) # here sqrt(dte) from A and from the SNR expression cancel 
  # as a result of normalizing A to have rel. error at t=0 rerr0
}
rError <- function(model, x, rerr0=3e-2) rerr0 * exp(-.5*model@beamLam * x) 
.ggplot_Sgl <- function(s, npts=300){
  ggplot(s[seq(1,nrow(s), length.out=npts),]) + 
    geom_line(aes(Time, XSgl), lwd=.2) + 
    geom_point(aes(Time, Sgl), col="red", size=.5) +
    thm +labs(y="Signal")
}

## SIGNAL ####
mod = CModel(); mod@beamLam <- 0#mod@decohLam
smpl = CuSampling()
taud = -1/mod@decohLam
taub = -1/mod@beamLam
tau = taud/(1+taud/taub)
Ttot=3*tau

s = simSample(smpl, mod, Ttot)
.ggplot_Sgl(s, 502) 

## COMPACTION TIME VS W SE ####
vwt = varwT(s, mod)
Gtot = g(Ttot)

SDe = smpl@rerror*mod@Num0

dte = 2e-4
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
cmpts <- seq(1, .2, -.1)
test <- function(cmpt, mod, smpl){
  smpl@CMPT <- cmpt
  mod@beamLam <- mod@beamLam*as.numeric(cmpt)
  
  w = mod@wFreq; dt = cmpt*pi/w
  
  taud = -1/mod@decohLam
  taub = -1/mod@beamLam
  tau = taud/(1+taud/taub)
  
  
  N0 <- mod@Num0; P <- mod@Pol; w0 <- mod@wFreq; phi <- mod@Phase
  lam.decoh <- mod@decohLam; lam.beam <- mod@beamLam
  
  f = Sgl ~ N0 * exp(lb*Time) * (1 + P*exp(ld*Time)*sin(w*Time + phi))
  guess = llply(c("w" = w0, "lb"= lam.beam, "ld" = lam.decoh), function(x) rnorm(1, x, abs(x)*1e-4))
  
  s = simSample(smpl, mod, 3*tau)
  nezc <- as.numeric(s[,.(Num=length(Time)), by=Node][5,Num])
  x0 = .5*(1+ sin(w*dt)/w/dt)
  G=g(s[length(Time),Time])
  
  nls(f, s, guess) %>% summary %>% coef -> stats
  
  
  list(Fit=stats, Smpl=c(
    Size=nrow(s), PhysT=s[length(Time),Time]-s[1,Time], Spread=sqrt(varwT(s, mod)),
    X0=x0, G=G, Nezc=nezc, FItot=x0*G*nezc
  ))

  
  # c(StatT=sqrt(varwT(s, mod)), Xtot=nezc*g(s[length(Time),Time])*.5*(1+ sin(w*dt)/w))
}

mod = CModel(); mod@beamLam <- mod@decohLam
smpl <- CmSampling(sglFreqGuess=rnorm(1,mod@wFreq,1e-4))

registerDoParallel(detectCores()-1)
llply(cmpts, function(x, m, s) test(x, m, s), mod, smpl,
      .parallel = TRUE, .paropts = list(.export=lsf.str(), .packages="dplyr")) -> stats
stopImplicitCluster()

stats.smpl <- ldply(stats, function(e) e$Smpl, .id="CMPT")
stats.fit <- ldply(stats, function(e) e$Fit["w",1:2], .id="CMPT")
stats.smpl<-mutate(stats.smpl, Denom = Spread*sqrt(FItot))
stats.smpl<-mutate(stats.smpl, 
                   Denom=Denom/Denom[1], 
                   FItot=FItot/FItot[1], 
                   Spread=Spread/Spread[1], 
                   Size=Size/Size[1],
                   X0=X0/X0[1],
                   Nezc=Nezc/Nezc[1]
                   )
stats.smpl<-melt(stats.smpl, id.vars = "CMPT", variable.name="parameter")

filter(stats.smpl, parameter %in% c("Size", "Spread", "FItot","Denom")) %>%
  ggplot(aes(CMPT, value, col=parameter)) + 
  geom_point() + 
  geom_line(aes(as.numeric(CMPT), value)) +
  thm

ggplot(stats.fit, aes(CMPT, `Std. Error`)) + geom_point() + thm + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

stats <- setDT(stats)
stats[,`:=`(Num=length(Time), Last=Time[length(Time)]), by=CMPT]

ggplot(stats[seq(1,nrow(stats), length.out = 2620),]) + geom_point(aes(Time, XSgl, col=CMPT),size=.2) +thm

#stats.w <- ldply(stats, function(e) e["w",1:2], .id="CMPT")
#ggplot(stats.w) + geom_pointrange(aes(CMPT, Estimate, ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`))

## comparison of modulated vs uniform samplings
smpl0 <- CmSampling(Freq=50, CMPT=.4)
smpl00 <- CmSampling(Freq=50, CMPT=1)
simSample(smpl00, mod, c(2,4.21))->s00
simSample(smpl0, mod, c(2,4.21))->s0
ggplot(s0) + geom_line(aes(Time, XSgl), data=s00, size=.2) + 
  geom_point(aes(Time, Sgl), data=s00, size=.5, col="darkgreen") +
  geom_point(aes(Time, Sgl), data=s0, size=.5, col="red") +
  geom_vline(aes(xintercept=Node), data=s0, size=.15, col="blue") +
  thm


## DURATION ####
test <- function(n, f, guess){
  s = simSample(smpl, mod, c(0, n)*tau)
  (nls(f, s, guess) %>% summary %>% coef)["w", 1:2]
}

mod <- CModel(); smpl <- CmSampling()

taud = -1/mod@decohLam
taub = -1/mod@beamLam
tau = taud/(1+taud/taub)

N0 <- mod@Num0; P <- mod@Pol; w0 <- mod@wFreq; phi <- mod@Phase
lam.decoh <- mod@decohLam; lam.beam <- mod@beamLam

f = Sgl ~ N0 * exp(lb*Time) * (1 + P*exp(ld*Time)*sin(w*Time + phi))
guess = llply(c("w" = w0, "lb"= lam.beam, "ld" = lam.decoh), function(x) rnorm(1, x, abs(x)*1e-4))

ntau=seq(1:6); names(ntau) <- ntau
registerDoParallel(detectCores()-1)
ldply(ntau, function(n,func,ini) test(n, func, ini), f, guess, 
      .parallel = TRUE,
      .paropts = list(.export=lsf.str(), .packages="dplyr"),
      .id="NTau"
      ) -> stats
stopImplicitCluster()

ggplot(stats%>%mutate(`Std. Error`=`Std. Error`/`Std. Error`[6]), aes(NTau, `Std. Error`)) +
  geom_point() + thm
