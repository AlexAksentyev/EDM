library(plyr)
library(ggplot2); library(gridExtra)
library(grid)

rm(list=ls(all=TRUE))

source("./RScripts/RCBunch.R")
source("./RScripts/RCSignal.R")

## COMPUTATIONS ####
bl <- replicate(8, RCBunch$new(Npart=1, SDdy=2e-2, SDphi=.5e-2))
sl <- llply(bl, function(b) RCSignal$new(b, seq(0,10, .25/b$Synch["wFreq"])))
names(sl) <- as.character(1:length(sl))

df = ldply(sl, function(s) s$Signal, .id = "Ptcl")
s = sl[[1]]+sl[[2]]+sl[[3]]+sl[[4]]+sl[[5]]+sl[[6]]+sl[[7]]+sl[[8]]

ggplot(df, aes(Time, Val)) + geom_line(aes(col=Ptcl)) + geom_line(data=s$Signal) + 
  theme_bw() +theme(legend.position="top")

rm(bl, sl, s, df)

b1 <- RCBunch$new(Npart=1e4)
Tstt=0; Ttot=2000; dt = .5/b1$Synch["wFreq"] # pi/w0 to satisfy the Nyquist condition
stime <- seq(Tstt, Ttot, length.out = 1e4)

s1 <- RCSignal$new(b1, seq(Tstt, Ttot, length.out=1e4))

s1$fit()
dttol=1e-6
s1$findPts(what="Node", w.guess = coef(s1$Model)[2], tol=dttol)
s1$specPts %>% ddply("Which", function(h){
  x = c(h$Time[1], h$Time[1:(nrow(h)-1)])
  mutate(h, DT = Time-x)
}) -> spts

fitstat <- coef(summary(s1$Model))[,1:2]
fitstat[1,1] <- -1/fitstat[1,1]; fitstat[1,2] <- fitstat[1,1]^2*fitstat[1,2]
rownames(fitstat) <- c("tau","w")
fitstat <- formatC(fitstat, 3,format="e")

## PLOTS ##
## particle distributions ####
.gghist_plot <- function(df, name){
  ggplot(df, aes_string(name)) +
    geom_histogram(aes(y=..density..), fill="white", color="black") +
    geom_density() + 
    # stat_function(fun=dnorm, 
    #               args=list(mean=attr(df[,name],"Synch"), 
    #                         sd=attr(df[,name],"SD")), 
    #               col="red", lty=2) +
    theme_bw()
}
.gghist_plot(b1$EnsPS, "wFreq") + labs(x=expression(omega)) -> whist
.gghist_plot(b1$EnsPS, "Phi") + labs(x=expression(phi))-> phist
grid.arrange(whist,phist)

## optimization function ####
fn <- function(x) colSums( sin(b1$Phase(x)) )
dx = pi/b1$Synch["wFreq"]
x = seq(spts[5,"Time"]-dx, spts[5,"Time"]+dx, length.out=250)
df = data.frame(Time=x, Sgl=fn(x), Tgt=fn(x)^2/1e3) %>% melt(id.vars="Time")

ggplot(df, aes(Time, value)) + geom_line(aes(linetype=variable)) +
  geom_hline(yintercept = 0, col="red") + theme_bw() + labs(y="Variable") +
  scale_linetype_manual(breaks=c("Sgl","Tgt"),values=c(1,2), labels=c("Signal","Target/1e3")) +
  theme(legend.position="top")

## freq creep ####
filter(spts,N>1, Which=="Optim") %>% mutate(w = pi/DT, SEw = pi*sqrt(2)*dttol/DT^2) %>% filter(w<6000, w>2)%>%
  ggplot(aes(Time, w)) + geom_linerange(aes(ymin=w-SEw, ymax=w+SEw), size=.1) + geom_hline(yintercept=b1$Synch["wFreq"], col="red") +
  # geom_hline(yintercept=coef(s1$Model)[2]) + 
  theme_bw() + labs(y=expression(omega(t)))

## signal ####
ggplot(s1$Signal, aes(Time, Val)) + geom_line(col="red",lwd=.05) + 
    theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
  theme(legend.position="top")+
  geom_point(aes(col=Which), size=1, data=s1$specPts, show.legend = TRUE) +
  scale_color_manual(values = c("black","blue")) +
  annotation_custom(tableGrob(fitstat),
                    xmin=1000,ymin=-1000, ymax=-650) -> sglplot

t0 = .7*Ttot; dt0=25
ggplot(filter(s1$Signal, Time>t0&Time<t0+dt0), aes(Time,Val)) + geom_line(col="red") + theme_bw() +
  scale_color_manual(values=c("black","blue"))+
  geom_point(aes(Time, Val, col=Which), size=1, data=filter(s1$specPts, Time>t0&Time<t0+dt0))+labs(x="",y="") -> sglslc

vp <- viewport(width = .4, height = .3, x = .6, y = .5, just = c("right","center"))
print(sglplot); print(sglslc,vp=vp)

## spectrum ####
s1$Spectrum()->fps
