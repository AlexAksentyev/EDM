library(nlreg)
library(ggplot2)

rm(list=ls(all=TRUE))

lblfnt = 20
thm = theme_bw() + theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
                         legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

source("./RScripts/CModel.R")
source("./RScripts/CSampling.R")

.ggplot_Sgl <- function(s, npts=300, variable="Sgl"){
  require(ggplot2)
  ggplot(s[seq(1,nrow(s), length.out=npts),]) + 
    geom_line(aes_string("Time", paste0("X",variable)), lwd=.2) + 
    geom_point(aes_string("Time", variable), col="red", size=.5) +
    thm +labs(y="Signal")
}
.SNR <- function(x, taub, taud) 33* exp(-x/taud * (1+.5*taud/taub))

L = CModel(Phase=-pi/2); L@beamLam <- log(.25)/1000
R = CModel(Phase=+pi/2); R@beamLam <- log(.25)/1000

smpl = CuSampling(Freq=375)
# smpl = CmSampling(sglFreqGuess=rnorm(1,3,1e-6), CMPT=.2, Freq=375/.2)

Ttot=1000
simSample(smpl, L, Ttot, grow=TRUE) -> Ls
simSample(smpl, R, Ttot, grow=TRUE) -> Rs
df = rbind(Ls%>%mutate(Detector="Left"), Rs%>%mutate(Detector="Right"))
ggplot(df[seq(1,nrow(df),length.out=500),]) + facet_grid(Detector~.) + thm +
  #geom_line(aes(Time, XSgl, col=Detector)) +
  # geom_point(aes(Time, Sgl), col="black",size=.2) +
  #labs(y="Counting rate, #", x="Time, sec")
   geom_line(aes(Time, (Sgl-XSgl)/XSgl, col=Detector), lwd=.2) +
   labs(y=expression(frac(tilde(N)[0](t)-N[0](t),N[0](t))), x="Time, sec")

Alr = merge(Ls[,FIDrvt:=NULL], Rs[,FIDrvt:=NULL], by="Time", suffixes = c(".L",".R")); rm(Ls, Rs)
Alr[,`:=`(XVal=(XSgl.L-XSgl.R)/(XSgl.L+XSgl.R), Val=(Sgl.L-Sgl.R)/(Sgl.L+Sgl.R))]
ggplot(Alr[seq(1, nrow(Alr), length.out = 433),]) + thm + 
  # geom_line(aes(Time, Val-XVal), lwd=.2) + labs(y="Error") #plot of residuals
   geom_line(aes(Time, XVal), lwd=.1, col="black") + geom_point(aes(Time, Val), size=.2, col="red") + labs(y="Asymmetry", x="Time, sec")# plot of signal


f = Val ~ n0*exp(lam*Time)*sin(w*Time + p)
guess = list(n0=.4, lam=-1e-4, w=3, p=-pi/2)
nlreg(f, weights=~1/(n0*exp(lam*Time)), data=Alr, start=guess) -> mod3l

mod3l%>%summary -> mod3l.s
print(mod3l.s)

