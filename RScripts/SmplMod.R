source("./RScripts/classes.R")
source("./RScripts/definitions.R")

mod = CModel()
smpl = CuSampling()
sew = 1e-2

Ttot = -1/mod@decohLam * 2
dtmeas = 5e-3
Nnd = mod@wFreq*Ttot/pi

a = mod@wFreq^2/12
err = smpl@rerror*mod@Pol*mod@Num0
df = simSample(smpl, mod, Ttot); .varWT(df) -> v
lam = mod@decohLam * pi/mod@wFreq
const = as.numeric((err/sew)^2/v["VarWT"]) * dtmeas * (exp(lam)-1)/(exp(lam*Nnd)-1)

f <- function(dt) a*dt^4 -dt^2 + const
x = seq(0,1.2,.1); y = f(x); plot(x,y, type="l"); abline(h=0, col="red")

uniroot(f, c(0,pi/(2*mod@wFreq)))$root -> cmptime
cmptime/(pi/mod@wFreq) * 100 # percents of compaction time per node time
