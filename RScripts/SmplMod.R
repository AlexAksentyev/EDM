source("./RScripts/classes.R")
source("./RScripts/definitions.R")

library(ggplot2)

mod = CModel()
smpl = CuSampling()
sew = 1e-3

Ttot = -1/mod@decohLam * 2.4
dtmeas = 5e-3
Nnd = mod@wFreq*Ttot/pi

err = smpl@rerror*mod@Pol*mod@Num0
df = simSample(smpl, mod, Ttot); .varWT(df) -> v
lam = mod@decohLam * pi/mod@wFreq
const = 2*as.numeric((err/sew)^2/v["VarWT"]) * dtmeas * (exp(lam)-1)/(exp(lam*Nnd)-1)

f <- function(dt) dt^2 + dt*sin(mod@wFreq*dt)/mod@wFreq - const
x = seq(-.1,1.5,.1); y = f(x); plot(x,y, type="l"); abline(h=0, col="red")

uniroot(f, c(0,2*pi/(2*mod@wFreq)))$root -> cmptime
cmptime/(pi/mod@wFreq) * 100 # percents of compaction time per node time

## plots
n=5 #number of seconds
ptime = smpl@Freq*n
nds = c(0, seq(n)*pi/mod@wFreq)
df[seq(1,ptime,length.out=250),]%>%ggplot(aes(Time, Sgl)) + geom_line() + theme_bw() +
  geom_vline(xintercept = -.5*cmptime + nds, col="blue") +
  geom_vline(xintercept = .5*cmptime + nds, col="red") +
  geom_point(aes(x,y), data.frame(x=nds, y=mod@Num0))
