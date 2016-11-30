library(dplyr); library(plyr)
fancy <- function(x) formatC(x, 4, format = "e")

## I should show that my formula for var omega is the correct one
## define the measurement model
P = .8; N0 = 6730
errS = 1e-1*N0*P #absolute measurement error

phi = pi/2 #initial signal phase
w = 2*pi*.005#2*pi*120e3*.74e-4 # estimated frequency
.dcs <- function(x) N0 * (1 + P*sin(w*x + phi)) # structural model/expectation

w.guess = rnorm(1, w, 0*w) # if we do modulated sampling, we must guess the signal frequency; this line models the guessing

Ttot = 1000 #total measurement time
fs = w/(2*pi)*sqrt(8) #sampling frequency
t = seq(0, Ttot, by = 1/fs);#uniform sampling time 
t2 = seq(0, Ttot, by = .33/fs); plot(.dcs(t2)~t2, type="o"); points(.dcs(t)~t, col="blue")
t2 <- t2[.dcs(t2) > N0*(1 - P*.5) & .dcs(t2) < N0*(1 + P*.5)]; points(.dcs(t2)~t2, col="red"); 
abline(h=N0*(1 + c(-.5,0,.5)*P), col="red", lty=c(2,3,2))
x = seq(0,Ttot,length.out=1000); lines(.dcs(x)~x, lty=4)
min(length(t), length(t2))->l
t <- t[1:l]; t2 = t2[1:l]
t <- list("uni" = t, "mod" = t2)
rm(t2)

Es = llply(t, function(x) {.dcs(x) -> e; names(e) <- x; e}) # expected values for the two samplings

Ntrial = 5 #number of trials
matrix(rnorm(Ntrial*l, sd=errS), ncol=Ntrial) -> s
llply(Es, function(e){e+s -> x; rownames(x) <- names(e); x})->s

## example of one trial
trl = llply(s, function(e) data.frame(t = as.numeric(rownames(e)), m = e[,3]))
ci = .95; z = abs(qnorm((1-ci)/2))
llply(trl, function(e){
  nls(m~ n*(1 + p*sin(frq*t + phi)), data=e, start=list(n = 5000, p = 1, frq=w.guess)) -> m3
  frq = coef(m3)["frq"]; sef = sqrt(vcov(m3)["frq","frq"])
  confint(m3) %>% cbind("SE" = sqrt(diag(vcov(m3))), "True" = c(N0, P, w)) %>% fancy
})



######
# adply(s, 2, function(trl){
#   nls(trl~ n*(1 + p*sin(frq*t)), start=list(n = 5000, p = 1, frq=rnorm(1, w, .1*w)))
# })