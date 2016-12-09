library(dplyr)

#### functions ####
fancy <- function(x) formatC(x, 4, format = "e")
.se <- function(x) sd(x)/sqrt(length(x))
.extract.stats <- function(m){
  .est = coef(m); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}

# 
# .dcs <- function(x, w = w0) N0 * (1 + P*exp(lam.decoh*x)*sin(w*x + phi)) # structural model/expectation
# .ddcs <- function(x, w = w0) N0*P*exp(lam.decoh*x)*cos(w*x + phi) # signal derivative

.compVarF <- function(df, err = errS){
  ftr <- sum(df$Drvt^2)
  mutate(df, Wt = Drvt^2/ftr, WtT = Time*Wt, WtTT = Time^2*Wt)->df
  (sum(df$WtTT) - sum(df$WtT)^2)*ftr -> denom
  err/sqrt(denom)
}

# .sample <- function(Nprd, Ntrl = 1){
#   Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
#   assign("Ttot", Ttot, envir = .GlobalEnv)
#   t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
#   
#   ## modulated sampling if the measurements are taken uniformly in time
#   # t2 = seq(0, Ttot, by = comptn/fs)
#   # x = .dcs(t2, w=w.g) #w.g reflects the fact that we don't know exactly the frequency we estimate;
#   #         # so, the points that pass the `goodness` test below pass it for the guessed signal,
#   #         #but not for the actual one
#   # DeltaS = P*comptn*exp(lam.decoh*t2)
#   # t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS)] #pick the good points
#   # rm(x) 
#   # min(length(t1), length(t2))->l
#   # t1 <- t1[1:l]; t2 <- t2[1:l]
#   ##
#   ## modulated sampling if the measurements are aggregated about the signal zero-crossing
#   ## the zero-crossing times are tn = (n*pi - phi)/w0 (but in our case we guess, so, w0 -> w.g)
#   ## the local sampling range should still correspond to the compaction factored signal range, 
#   ## this determines Dt; Dt'll have to change with time, if we are to keep the modulation efficiency gain
#   ## !!!! this will have to be done analytically, but for now I'll just 
#   Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)] #time range about the z-crossings
#   tn = ((0:(2*Nprd))*pi - phi)/w.g # guessed z-crossing times
#   t2 = laply(tn, function(ti) ti+Dt) %>% c 
#   x = .dcs(t2);  DeltaS = P*comptn*exp(lam.decoh*t2) # information condition
#   t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0] #picking the good points
#   rm(x) 
#   min(length(t1), length(t2))->l
#   t1 <- t1[1:l]; t2 <- t2[order(t2)][1:l]
#   ##
#   
#   df1 = data.frame("Type" = "uniform","Time" = t1, "XSgl" = .dcs(t1), "Drvt" = .ddcs(t1))
#   df2 = data.frame("Type" = "modulated","Time" = t2, "XSgl" = .dcs(t2), "Drvt" = .ddcs(t2))
#   rbind(df1,df2)[rep(seq(2*l),Ntrl),] %>% arrange(Type) %>%
#     mutate(Trl = rep(seq(Ntrl), each=l)%>%rep(2), Sgl = XSgl + rnorm(Ntrl*l, sd=errS)%>%rep(2))
# }
.fit <- function(.sample, model, return.model = FALSE){
  w.g = rnorm(1, model@wFreq, .001)
  n.g = rnorm(1, model@Num0, .1*model@Num0)
  p.g = rnorm(1, model@Pol, .1*model@Pol)
  lam.decoh = model@decohLT
  phi = model@Phase
  
  cat("guessed frequency:", w.g, "\n\n")
  
  nls(Sgl ~ n*(1 + p*exp(lam.decoh*Time)*sin(frq*Time + phi)), data=.sample, start=list(n = n.g, p = p.g, frq = w.g)) -> m3
  
  if(return.model) return(m3) else return(.extract.stats(m3))
}

varW0_test <- function(.mod, .Time){
  # cat(paste("model frequency", .mod$wfreq, "\n"))
  
  simSample(stu, .mod, .Time) -> .spl
  # cat(paste("sample size", nrow(.spl), "\n"))
  
  .fit(.spl, .mod) -> .stats
  # .stats = NULL
  
  list("Stats" = .stats, "Sample" = .spl)
}

#only modulated sampling
# .msampleF <- function(Nprd = 5, fs = 5000, len = NA, comptn = .5){
#   Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
#   assign("Ttot",Ttot, envir = .GlobalEnv)
#   Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)] #time range about the z-crossings
#   tn = ((0:(2*Nprd))*pi - phi)/w.g # guessed z-crossing times
#   t2 = laply(tn, function(ti) ti+Dt) %>% c 
#   x = .dcs(t2);  DeltaS = P*comptn*exp(lam.decoh*t2) # information condition
#   t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0] #picking the good points
#   rm(x)
#   if(is.na(len)) len <- length(t2)
#   t2 <- t2[order(t2)][1:len]
#   
#   data.frame("Time" = t2, "XSgl" = .dcs(t2), "Drvt" = .ddcs(t2)) -> df
#   df %>%  mutate(Sgl = XSgl + rnorm(len, sd=errS))
# }
# #only uniform sampling
# .usampleF <- function(Ttot, fs = 5000){
#   t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
#   
#   data.frame("Time" = t1, "XSgl" = .dcs(t1), "Drvt" = .ddcs(t1)) -> df
#   df %>%  mutate(Sgl = XSgl + rnorm(length(t1), sd=errS))
# }

#### classes ####
model <- function(wfreq=3, phs=pi/36, N0=6730, Pol=.4, decoh.LT = log(.25)/1000){

  me <- list(
    wfreq = wfreq, phase = phs,
    Num0 = N0, Pol = Pol,
    lam.decoh = decoh.LT
  )
  
  class(me) <- append(class(me), "model")
  return(me)
  
}
msampling <- function(freq=5000, freq.guess = NA, compaction =.3){
  guess = ifelse(is.na(freq.guess), rnorm(1, 3, .003), freq.guess)
  me <- list(
    freq = freq, cmptn = compaction,
    sglWfreq.guess = guess
  )
  
  class(me) <- append(class(me), "msampling")
  return(me)
}
usampling <- function(freq=5000){
  me <- list(freq = freq)
  
  class(me) <- append(class(me), "usampling")
  return(me)
}
# errS = 3e-2*N0*P #absolute measurement error

#### methods ####
setWFreq <- function(whose, value) UseMethod("setWFreq", whose)
setWFreq <- function(whose, value) {whose$wfreq <- value; whose}
setDecoh <- function(whose, value) UseMethod("setDecoh", whose)
setDecoh <- function(whose, value) {whose$lam.decoh <- value; whose}
expectation <- function(whose, which.pts) UseMethod("expectation",whose)
expectation.model <- function(model, at.points){
  N0 <- model$Num0; P <- model$Pol; w <- model$wfreq; phi <- model$phase
  lam.decoh <- model$lam.decoh
  
  N0 * (1 + P*exp(lam.decoh*at.points)*sin(w*at.points + phi))
}
derivative <- function(whose, where) UseMethod("derivative", whose)
derivative.model <- function(model, at.points){
  N0 <- model$Num0; P <- model$Pol; 
  w <- model$wfreq; phi <- model$phase
  lam.decoh = model$lam.decoh
  
  N0*P*exp(lam.decoh*at.points)*cos(w*at.points + phi)
  
}

sample <- function(how, model, how.long, err = NA) UseMethod("sample", how)
sample.usampling <- function(how, model, how.long, err = NA){
  if(is.na(err)) err <- 3e-2 * model$Num0*model$Pol
  
  t1 = seq(0, how.long, by = 1/how$freq) #uniform sampling
  
  data.frame("Time" = t1, "XSgl" = expectation(model, t1), "Drvt" = derivative(model, t1)) %>%  
    mutate(Sgl = XSgl + rnorm(length(t1), sd=err))
}
sample.msampling <- function(how, model, how.long, err = NA){
  if(is.na(err)) err <- 3e-2 * model$Num0*model$Pol
  
  phi = model$phase; w0 = model$wfreq; lam.decoh = model$lam.decoh; 
  P = model$Pol; N0 = model$Num0
  fs = how$freq; wg = how$sglWfreq
  cptn = how$cmptn
  
  Nprd = round((how.long*w0 + phi)/(2*pi)); cat(paste("periods", Nprd, "\n"))
  Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)]
  tn = ((0:(2*Nprd))*pi - phi)/wg; cat(paste("last z-c", tn[length(tn)], "\n"))
  t2 = laply(tn, function(ti) ti+Dt) %>% c 
  
  x = expectation(model, t2);  DeltaS = P*cptn*exp(lam.decoh*t2) # information condition
  t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0]
  t2 <- t2[order(t2)]
  
  data.frame("Time" = t2, "XSgl" = expectation(model, t2), "Drvt" = derivative(model, t2)) %>%
    mutate(Sgl = XSgl + rnorm(length(t2), sd=err))
}