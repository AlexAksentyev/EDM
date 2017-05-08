source("./RScripts/RCBunch.R")
source("./RScripts/RCSignal.R")

b = RCBunch$new()
tau = 5; dt = .5*pi/b$Synch["wFreq"]
s = RCSignal$new(b, seq(0,tau,dt))

for(i in 1:floor(1000/tau)){
  b$flipFreq(i*tau+dt)
  s$sample(seq(i*tau+dt,(i+1)*tau,dt))
}

ggplot(s$Signal, aes(Time, Val, col=Smpl)) + geom_point() + scale_color_continuous(low="blue",high="red")
