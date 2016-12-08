library(dplyr); library(plyr)

#### functions ####
fancy <- function(x) formatC(x, 4, format = "e")
.se <- function(x) sd(x)/sqrt(length(x))

.dcs <- function(x, w = w0) N0 * (1 + P*exp(lam.decoh*x)*sin(w*x + phi)) # structural model/expectation
.ddcs <- function(x, w = w0) N0*P*exp(lam.decoh*x)*cos(w*x + phi) # signal derivative

.compVarF <- function(df, err = errS){
  ftr <- sum(df$Drvt^2)
  mutate(df, Wt = Drvt^2/ftr, WtT = Time*Wt, WtTT = Time^2*Wt)->df
  (sum(df$WtTT) - sum(df$WtT)^2)*ftr -> denom
  err/sqrt(denom)
}

.sample <- function(Nprd, Ntrl = 1){
  Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
  assign("Ttot", Ttot, envir = .GlobalEnv)
  t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
  
  ## modulated sampling if the measurements are taken uniformly in time
  # t2 = seq(0, Ttot, by = comptn/fs)
  # x = .dcs(t2, w=w.g) #w.g reflects the fact that we don't know exactly the frequency we estimate;
  #         # so, the points that pass the `goodness` test below pass it for the guessed signal,
  #         #but not for the actual one
  # DeltaS = P*comptn*exp(lam.decoh*t2)
  # t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS)] #pick the good points
  # rm(x) 
  # min(length(t1), length(t2))->l
  # t1 <- t1[1:l]; t2 <- t2[1:l]
  ##
  ## modulated sampling if the measurements are aggregated about the signal zero-crossing
  ## the zero-crossing times are tn = (n*pi - phi)/w0 (but in our case we guess, so, w0 -> w.g)
  ## the local sampling range should still correspond to the compaction factored signal range, 
  ## this determines Dt; Dt'll have to change with time, if we are to keep the modulation efficiency gain
  ## !!!! this will have to be done analytically, but for now I'll just 
  Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)] #time range about the z-crossings
  tn = ((0:(2*Nprd))*pi - phi)/w.g # guessed z-crossing times
  t2 = laply(tn, function(ti) ti+Dt) %>% c 
  x = .dcs(t2);  DeltaS = P*comptn*exp(lam.decoh*t2) # information condition
  t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0] #picking the good points
  rm(x) 
  min(length(t1), length(t2))->l
  t1 <- t1[1:l]; t2 <- t2[order(t2)][1:l]
  ##
  
  df1 = data.frame("Type" = "uniform","Time" = t1, "XSgl" = .dcs(t1), "Drvt" = .ddcs(t1))
  df2 = data.frame("Type" = "modulated","Time" = t2, "XSgl" = .dcs(t2), "Drvt" = .ddcs(t2))
  rbind(df1,df2)[rep(seq(2*l),Ntrl),] %>% arrange(Type) %>%
    mutate(Trl = rep(seq(Ntrl), each=l)%>%rep(2), Sgl = XSgl + rnorm(Ntrl*l, sd=errS)%>%rep(2))
}
.fit <- function(trl, model = FALSE){
  cat(paste("freq guess", w.g, "\n"))
  nls(Sgl~ n*(1 + p*exp(lam.decoh*Time)*sin(frq*Time + phi)), data=trl, start=list(n = 5000, p = 1, frq=w.g)) -> m3
  if(model) return(m3)
  
  .extract.stats(m3)
}
.extract.stats <- function(m){
  .est = coef(m); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}
#only modulated sampling
.msampleF <- function(Nprd = 5, fs = 5000, len = NA, comptn = .5){
  Ttot = (Nprd*2*pi-phi)/w0 #total measurement time
  assign("Ttot",Ttot, envir = .GlobalEnv)
  Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)] #time range about the z-crossings
  tn = ((0:(2*Nprd))*pi - phi)/w.g # guessed z-crossing times
  t2 = laply(tn, function(ti) ti+Dt) %>% c 
  x = .dcs(t2);  DeltaS = P*comptn*exp(lam.decoh*t2) # information condition
  t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0] #picking the good points
  rm(x)
  if(is.na(len)) len <- length(t2)
  t2 <- t2[order(t2)][1:len]
  
  data.frame("Time" = t2, "XSgl" = .dcs(t2), "Drvt" = .ddcs(t2)) -> df
  df %>%  mutate(Sgl = XSgl + rnorm(len, sd=errS))
}
#only uniform sampling
.usampleF <- function(Ttot, fs = 5000){
  t1 = seq(0, Ttot, by = 1/fs) #uniform sampling
  
  data.frame("Time" = t1, "XSgl" = .dcs(t1), "Drvt" = .ddcs(t1)) -> df
  df %>%  mutate(Sgl = XSgl + rnorm(length(t1), sd=errS))
}

#### parameters ####
w0=3.5; phi=pi/36; N0=6730; P=.4 ; lam.decoh = log(.25)/1000# signal; /25 b/c I want to model 1000 secs by 25 secs (12 periods)
fs=5000; comptn=.3; w.g=rnorm(1,w0,.001*w0) ## sampling; we guess the true frequency with 1% precision
errS = 3e-2*N0*P #absolute measurement error