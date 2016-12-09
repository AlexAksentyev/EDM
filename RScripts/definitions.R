library(dplyr)

#### utility functions ####

.se <- function(x) sd(x)/sqrt(length(x))
.extract.stats <- function(m){
  .est = coef(m); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}

.compVarF <- function(df, err = 3e-2){
  err <- err*df$XSgl[1]
  ftr <- sum(df$fishInfo^2)
  mutate(df, Wt = fishInfo^2/ftr, WtT = Time*Wt, WtTT = Time^2*Wt)->df
  (sum(df$WtTT) - sum(df$WtT)^2)*ftr -> denom
  err/sqrt(denom)
}

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

#### tests ####
#### 5) SE(estimate(freq)) vs freq ####
varW0_test <- function(model, sampling, duration = 100, wfreqs = c(.01,.1,.3,1,3,10)){
  require(parallel); require(doParallel); n.cores = detectCores()-1
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  rtns <- c(".extract.stats", ".fit", ".compVarF", "simSample","fisherInfo","expectation")
  clusterExport(clus, rtns)
  
  names(wfreqs) <- wfreqs
  llply(wfreqs, function(w) setValue(model, c("wFreq" = w))) %>%
    llply(
    function(.mod, .splstg, .Time){
      simSample(.splstg, .mod, .Time) -> .spl; .fit(.spl, .mod) %>%mutate(SEAN.frq=.compVarF(.spl))-> .stats; 
      list("Stats" = .stats, "Sample" = .spl)
    }, .splstg = sampling, .Time = duration, 
    .parallel=TRUE, 
    .paropts = list(.export=c(rtns), .packages=c("dplyr"))
  ) -> dat
  
  dat
  
}


#### S3 classes ####
# model <- function(wfreq=3, phs=pi/36, N0=6730, Pol=.4, decoh.LT = log(.25)/1000){
# 
#   me <- list(
#     wfreq = wfreq, phase = phs,
#     Num0 = N0, Pol = Pol,
#     lam.decoh = decoh.LT
#   )
#   
#   class(me) <- append(class(me), "model")
#   return(me)
#   
# }
# msampling <- function(freq=5000, freq.guess = NA, compaction =.3){
#   guess = ifelse(is.na(freq.guess), rnorm(1, 3, .003), freq.guess)
#   me <- list(
#     freq = freq, cmptn = compaction,
#     sglWfreq.guess = guess
#   )
#   
#   class(me) <- append(class(me), "msampling")
#   return(me)
# }
# usampling <- function(freq=5000){
#   me <- list(freq = freq)
#   
#   class(me) <- append(class(me), "usampling")
#   return(me)
# }
# # errS = 3e-2*N0*P #absolute measurement error
# 
# #### methods ####
# setWFreq <- function(whose, value) UseMethod("setWFreq", whose)
# setWFreq <- function(whose, value) {whose$wfreq <- value; whose}
# setDecoh <- function(whose, value) UseMethod("setDecoh", whose)
# setDecoh <- function(whose, value) {whose$lam.decoh <- value; whose}
# expectation <- function(whose, which.pts) UseMethod("expectation",whose)
# expectation.model <- function(model, at.points){
#   N0 <- model$Num0; P <- model$Pol; w <- model$wfreq; phi <- model$phase
#   lam.decoh <- model$lam.decoh
#   
#   N0 * (1 + P*exp(lam.decoh*at.points)*sin(w*at.points + phi))
# }
# derivative <- function(whose, where) UseMethod("derivative", whose)
# derivative.model <- function(model, at.points){
#   N0 <- model$Num0; P <- model$Pol; 
#   w <- model$wfreq; phi <- model$phase
#   lam.decoh = model$lam.decoh
#   
#   N0*P*exp(lam.decoh*at.points)*cos(w*at.points + phi)
#   
# }
# 
# sample <- function(how, model, how.long, err = NA) UseMethod("sample", how)
# sample.usampling <- function(how, model, how.long, err = NA){
#   if(is.na(err)) err <- 3e-2 * model$Num0*model$Pol
#   
#   t1 = seq(0, how.long, by = 1/how$freq) #uniform sampling
#   
#   data.frame("Time" = t1, "XSgl" = expectation(model, t1), "Drvt" = derivative(model, t1)) %>%  
#     mutate(Sgl = XSgl + rnorm(length(t1), sd=err))
# }
# sample.msampling <- function(how, model, how.long, err = NA){
#   if(is.na(err)) err <- 3e-2 * model$Num0*model$Pol
#   
#   phi = model$phase; w0 = model$wfreq; lam.decoh = model$lam.decoh; 
#   P = model$Pol; N0 = model$Num0
#   fs = how$freq; wg = how$sglWfreq
#   cptn = how$cmptn
#   
#   Nprd = round((how.long*w0 + phi)/(2*pi)); cat(paste("periods", Nprd, "\n"))
#   Dt = c(seq(-3,3,2/fs), seq(-1.5,1.5,1/fs)); Dt <- Dt[order(Dt)]
#   tn = ((0:(2*Nprd))*pi - phi)/wg; cat(paste("last z-c", tn[length(tn)], "\n"))
#   t2 = laply(tn, function(ti) ti+Dt) %>% c 
#   
#   x = expectation(model, t2);  DeltaS = P*cptn*exp(lam.decoh*t2) # information condition
#   t2 <- t2[x > N0*(1 - DeltaS) & x < N0*(1 + DeltaS) & t2 >= 0]
#   t2 <- t2[order(t2)]
#   
#   data.frame("Time" = t2, "XSgl" = expectation(model, t2), "Drvt" = derivative(model, t2)) %>%
#     mutate(Sgl = XSgl + rnorm(length(t2), sd=err))
# }