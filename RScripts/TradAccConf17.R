source("./RScripts/classes.R")
source("./RScripts/definitions.R")

library(ggplot2)
library(reshape2)

lblfnt=20
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

mod = CModel()
smpl = CuSampling()

## LOOKING FOR TOTAL MEASUREMENT TIME #### 
## THIS TIME, VIA AN INFORMATIVITY CRITERION
Ttot = 1800
g <- function(x, dlt, limit = FALSE){
  lam = -1/dlt; denom = exp(lam/mod@wFreq * pi)-1; 
  
  if(limit) return (- 1/denom)
  else return ((exp(lam*x)-1)/denom)
} 

x = 0:Ttot; LTs = c("1000" = 1000, "721" = 721, "500" = 500, "250" = 250)
g(0, LTs, TRUE) -> inflims
ldply(0:Ttot, function(x) c("Time" = x, g(x, LTs))) %>% melt(id.vars="Time", variable.name="dLT", value.name = "g") -> dat

ggplot(dat) + geom_line(aes(Time, g, col=dLT)) + thm +
  geom_hline(yintercept = inflims, lty=2) +
  scale_color_discrete(guide = guide_legend(title=expression(tau~" "))) +
  scale_x_continuous(name="Measurement time (s)") +
  scale_y_continuous(name="")

## MODULATION ####
err = smpl@rerror*mod@Pol*mod@Num0
lam = mod@decohLam * pi/mod@wFreq

sew = c(c(1,2,5)%o%10^(-4:-2))
dtmeas = 5e-3

## VarWT for the different experiment durations
library(doParallel)
makeCluster(detectCores()) -> clus; registerDoParallel(clus)
rtns <- lsf.str(envir=.GlobalEnv, all=TRUE)
clusterExport(clus, rtns)

part = c(.7, 1.2, 2.3, 3.0)
Ttot = -1/mod@decohLam*part; names(Ttot) <- as.character(part)
Nnd = floor(mod@wFreq*Ttot/pi)
laply(Ttot, function(x) .varWT(simSample(smpl, mod, x))["VarWT"], 
      .parallel = TRUE, .paropts = list(.export=c("smpl","mod"), .packages="dplyr"))%>%unlist->v
stopCluster(clus)
names(v) <- as.character(part)

## the relevant constant values
const = 2*outer((err/sew)^2,v,"/") * dtmeas * (exp(lam)-1)/(exp(lam*Nnd)-1); rownames(const)<-as.character(sew)

## searching for the solutions
adply(const, c(1,2), function(a) 
  tryCatch(uniroot(function(dt) dt + sin(mod@wFreq*dt)/mod@wFreq - a, c(0,pi/mod@wFreq))$root, 
           error=function(e) pi/mod@wFreq)
) %>% `names<-`(c("SEW", "Part", "CMPT"))-> cmptime

library(mosaic)
mutate(cmptime, fCMPT = derivedFactor(
  ">50%" = CMPT/pi*mod@wFreq>.5,
  ">40%" = CMPT/pi*mod@wFreq>.4,
  ">30%" = CMPT/pi*mod@wFreq>.3,
  ">20%" = CMPT/pi*mod@wFreq>.2,
  ">10%" = CMPT/pi*mod@wFreq>.1,
  .default = "<10%",
  .method = "first"
)) -> cmptime

## plots 
ggplot(cmptime) + geom_tile(aes(SEW,Part,fill=fCMPT)) + 
  labs(x=expression(sigma[hat(omega)]), y=expression(""%*%tau[d])) + thm +
  scale_fill_discrete(guide=guide_legend(title="")) +
  scale_x_discrete(breaks=sew[seq(1,length(sew),length.out = 5)], labels=.fancy_scientific)

## STATISTICAL PRECISION ####
msmpl = CmSampling(sglFreqGuess=rnorm(1,3,3e-4), CMPT=.2)
simSample(msmpl, mod, Ttot["2.3"]) -> s
.ggplot_XSmpl(s, mod) + #geom_point(aes(Node, mod@Num0, col="red"), show.legend = FALSE) +
  theme_bw()

# part = c(.7, 1.2, 2.3, 3)
part=seq(.5,.9,.1)
Ttot = -1/mod@decohLam*part; names(Ttot) <- as.character(part)
expand.grid(Ttot=Ttot, CMPT = c(.05,.2,.3,1)*pi/mod@wFreq) %>% 
  adply(1, function(p, sampling, model){
    s=setValue(sampling, c("CMPT"=p[,"CMPT"]))
    simSample(s, model, p[,"Ttot"]) %>%
      .compAnaWSE(aerr=sampling@rerror*model@Num0)%>% `names<-`("SE")
    }, msmpl, mod
  ) -> x

mutate(x, fCMPT=as.factor(100*CMPT/pi*mod@wFreq%>%round(2))) %>% 
  ggplot(aes(as.numeric(names(Ttot)), SE*1e6, col=fCMPT)) + 
  geom_point() + geom_line() + thm +
  labs(x=expression(""%*%tau[d]), y=expression(sigma[hat(omega)]%*%10^6)) +
  scale_color_discrete(name="", labels=c("5%","20%","30%","100%"))


## COMBINED DECOHERENCE+SCATTERING LT ####
r = seq(.1,10,.25) #beam-LT/decoherence-LT
rTauSum = r/(1+r) # -(1/tau_b + 1/tau_d) / decoherence-LT
df=data.frame(R=r, rTS=rTauSum)
plot(df$R~df$rTS, type="l")

## RELATIVE ERROR OF THE COUTING RATE ####
p = 1e-2 # the fraction of the useful beam scatterings
nu = 1e6 # beam revolution frequency
N0b = 5e11 # number of particles in one fill
taud = -1/mod@decohLam #decoherence life-time
taub = taud * 2 # beam lifetime
lam = -1/taub-1/taud; tau = -1/lam
dtc = 1/nu # polarimetry measurement duration
dte = 2000*dtc # signal measurement duration

B = p*nu*N0b * 1e-10
N0 <- function(x) B*dtc*exp(-x/taub)
osc <- function(x, ampl=FALSE) {
  a = mod@Pol*exp(-x/taud); 
  if(ampl) return(a) 
  else return(a*sin(mod@wFreq*x+mod@Phase))
}
SEN0 <- function(x) sqrt(B/dte) * dtc * exp(-.5*x/taub)

df = data.frame(Time=c(.5,.7,1.2,2.3,3)*tau) %>% 
  mutate(rTime=Time/tau,
         No = N0(Time), Sgl = No*(1+osc(Time)), 
         SENo = SEN0(Time), rSENo=SENo/(No*(1+osc(Time,TRUE))), 
         SNR = No*osc(Time,TRUE)/SENo
  )
ggplot(df, aes(Time/tau, SNR)) + geom_line() + geom_point() + theme_bw() #+ geom_smooth(se=FALSE,col="red")

