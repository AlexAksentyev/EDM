library(dplyr); library(plyr)
library(ggplot2); library(gridExtra)
library(grid)

rm(list=ls(all=TRUE))

source("./RScripts/RCBunch.R")

## FUNCTIONS ####
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

## COMPUTATIONS ####
b1 <- RCBunch$new()

Tstt = 0; Ttot=2000; dt = .25/b1$Synch["wFreq"] # pi/w0 to satisfy the Nyquist condition
b1$project(seq(Tstt, Ttot, dt))
b1$fit()
b1$findPts("Envelope")
b1$Spectrum()

## PLOTS ####
.gghist_plot(b1$PS, "wFreq") + labs(x=expression(omega)) -> whist
.gghist_plot(b1$PS, "Phi") + labs(x=expression(phi))-> phist
grid.arrange(whist,phist)

ggplot(b1$Pproj, aes(Time, Val)) + geom_line(col="red",lwd=.05) + 
  theme_bw() + labs(y=expression(pi[bold(y)]*bold(P)))+
  theme(legend.position="top")+
  geom_point(aes(col=Which), size=1, data=b1$specPts, show.legend = TRUE) +
  scale_color_manual(values = c("black","blue")) +
  annotation_custom(tableGrob(formatC(coef(summary(b1$Model))[,1:2],3,format="e")), 
                    xmin=1000,ymin=-1000, ymax=-650) -> sglplot

t0 = 1900; dt0=25
ggplot(filter(b1$Pproj, Time>t0&Time<t0+dt0), aes(Time,Val)) + geom_line(col="red") + theme_bw() +
  scale_color_manual(values=c("black","blue"))+
  geom_point(aes(Time, Val, col=Which), size=1, data=filter(b1$specPts, Time>t0&Time<t0+dt0))+labs(x="",y="") -> sglslc

vp <- viewport(width = .4, height = .3, x = .6, y = .5, just = c("right","center"))
print(sglplot); print(sglslc,vp=vp)
