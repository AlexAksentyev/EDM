library(dplyr); library(plyr)
library(ggplot2); library(gridExtra)
library(grid)
library(reshape2)

rm(list=ls(all=TRUE))

source("./RScripts/RCBunch.R")

## COMPUTATIONS ####
b1 <- RCBunch$new(WDist="norm")

Tstt=7500; Ttot=12000; dt = .5/b1$Synch["wFreq"] # pi/w0 to satisfy the Nyquist condition
b1$project(seq(Tstt, Ttot, dt))

b1$fit()
b1$findPts(what="Node", w.guess = coef(b1$Model)[2], tol=1e-6) # w.guess = coef(b1$Model)[2]
b1$specPts %>% ddply("Which", function(h){
  x = c(h$Time[1], h$Time[1:(nrow(h)-1)])
  mutate(h, DT = Time-x)
}) -> spts

fitstat <- coef(summary(b1$Model))[,1:2]
fitstat[1,1] <- -1/fitstat[1,1]; fitstat[1,2] <- fitstat[1,1]^2*fitstat[1,2]
rownames(fitstat) <- c("tau","w")
fitstat <- formatC(fitstat, 3,format="e")

## PLOTS ##
## particle distributions ####
.gghist_plot <- function(df, name){
  ggplot(df, aes_string(name)) +
    geom_histogram(aes(y=..density..), fill="white", color="black") +
    geom_density() + 
    stat_function(fun=dnorm, 
                  args=list(mean=attr(df[,name],"Synch"), 
                            sd=attr(df[,name],"SD")), 
                  col="red", lty=2) +
    theme_bw()
}
.gghist_plot(b1$PS, "wFreq") + labs(x=expression(omega)) -> whist
.gghist_plot(b1$PS, "Phi") + labs(x=expression(phi))-> phist
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
filter(spts,N>1, Which=="Optim") %>% mutate(w = pi/DT) %>% filter(w<6000, w>2)%>%
  ggplot(aes(Time, w)) + geom_point(size=.3) + geom_hline(yintercept=b1$Synch["wFreq"], col="red") +
  geom_hline(yintercept=coef(b1$Model)[2]) + theme_bw() + labs(y=expression(omega(t)))

## signal ####
ggplot(b1$Pproj, aes(Time, Val)) + geom_line(col="red",lwd=.05) + 
    theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
  theme(legend.position="top")+
  geom_point(aes(col=Which), size=1, data=b1$specPts, show.legend = TRUE) +
  scale_color_manual(values = c("black","blue")) +
  annotation_custom(tableGrob(fitstat),
                    xmin=1000,ymin=-1000, ymax=-650) -> sglplot

t0 = 950; dt0=25
ggplot(filter(b1$Pproj, Time>t0&Time<t0+dt0), aes(Time,Val)) + geom_line(col="red") + theme_bw() +
  scale_color_manual(values=c("black","blue"))+
  geom_point(aes(Time, Val, col=Which), size=1, data=filter(b1$specPts, Time>t0&Time<t0+dt0))+labs(x="",y="") -> sglslc

vp <- viewport(width = .4, height = .3, x = .6, y = .5, just = c("right","center"))
print(sglplot); print(sglslc,vp=vp)

## spectrum ####
b1$Spectrum()->fps
