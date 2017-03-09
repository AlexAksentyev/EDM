source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

smpl = CuSampling(Freq=50)
L = CModel(Phase=-pi/2)
R = CModel(Phase=pi/2)

simSample(smpl, L, 500, grow=FALSE) -> Ls
simSample(smpl, R, 500, grow=TRUE) -> Rs
ggplot(Ls) + geom_line(aes(Time, XSgl), lwd=.2, col="blue") #+ geom_point(aes(Time, Sgl), size=.5)
