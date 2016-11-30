library(dplyr); library(plyr)

## I should show that my formula for var omega is the correct one

fancy <- function(x) formatC(x, 4, format = "e")

## define the measurement model
P = .8; N0 = 6730
errS = 1e-1*N0*P #absolute measurement error

phi = pi/2 #initial signal phase
w = 2*pi*.05#*120e3*.74e-4 # estimated frequency
w.guess = rnorm(1, w, .1*w)

fs = w.guess/2.33 #sampling frequency
t = seq(0, (pi-phi)/w.guess*Nprd, by = 1/fs) #sampling time
Nprd = 10
t <- list(
  "mod" = rep((pi-phi)/w.guess * 0:(Nprd-1), each = round(length(t)/Nprd)) + t[1:round(length(t)/Nprd)], 
  "uni" = t
)

Es = ldply(t, function(x) N0 * (1 + sin(w*x + phi))) #the structural model/expectation

Ntrial = 5 #number of trials
matrix(rnorm(Ntrial*length(t), sd=errS) + Es, ncol=Ntrial) -> s

## example of one trial
trl = s[,3]; plot(trl~t)
nls(trl~ n*(1 + p*sin(frq*t + phi)), start=list(n = 5000, p = 1, frq=w.guess)) -> m3
summary(m3); frq = coef(m3)["frq"]; sef = sqrt(vcov(m3)["frq","frq"])
ci = .95; z = abs(qnorm((1-ci)/2))
confint(m3) %>% cbind("SE" = sqrt(diag(vcov(m3))), "True" = c(N0, P, w)) %>%  print()


######
# adply(s, 2, function(trl){
#   nls(trl~ n*(1 + p*sin(frq*t)), start=list(n = 5000, p = 1, frq=rnorm(1, w, .1*w)))
# })