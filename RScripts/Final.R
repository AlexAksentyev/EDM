library(nlreg)

rm(list=ls(all=TRUE))

source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

.ggplot_Sgl <- function(s, npts=300, variable="Sgl"){
  require(ggplot2)
  ggplot(s[seq(1,nrow(s), length.out=npts),]) + 
    geom_line(aes_string("Time", paste0("X",variable)), lwd=.2) + 
    geom_point(aes_string("Time", variable), col="red", size=.5) +
    thm +labs(y="Signal")
}

smpl = CuSampling(Freq=500)
L = CModel(Phase=-pi/2); L@beamLam <- L@decohLam
R = CModel(Phase=pi/2); R@beamLam <- R@decohLam


Ttot=1082
simSample(smpl, L, Ttot, grow=TRUE) -> Ls
simSample(smpl, R, Ttot, grow=TRUE) -> Rs
df = rbind(Ls%>%mutate(Which="L"), Rs%>%mutate(Which="R")); rm(Ls, Rs)
ggplot(df[seq(1,nrow(df), length.out = 1200),]) + facet_grid(Which~.) + 
  geom_line(aes(Time, Sgl-XSgl, col=Which), lwd=.2) + theme(legend.position="top")

Alr = merge(Ls[,FIDrvt:=NULL], Rs[,FIDrvt:=NULL], by="Time", suffixes = c(".L",".R"))
Alr[,`:=`(XVal=(XSgl.L-XSgl.R)/(XSgl.L+XSgl.R), Val=(Sgl.L-Sgl.R)/(Sgl.L+Sgl.R))]
ggplot(Alr[seq(1, nrow(Alr), length.out = 4273),]) + #geom_line(aes(Time, Val-XVal), lwd=.2) #plot of residuals
  geom_line(aes(Time, XVal), lwd=.2, col="red") + geom_point(aes(Time, Val), size=.2)


f = Val ~ n0*exp(lam*Time)*sin(w*Time + p)
guess = list(n0=.4, lam=-1e-4, w=3, p=-pi/2)
nlreg(f, weights=~1/(n0*exp(lam*Time)), data=Alr, start=guess) -> mod3l

mod3l%>%summary -> mod3l.s
print(mod3l.s)
