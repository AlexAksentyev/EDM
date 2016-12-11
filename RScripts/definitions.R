library(dplyr)

#### utility functions ####
.fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
.se <- function(x) sd(x)/sqrt(length(x))
.extract.stats <- function(m){
  .est = coef(m); names(.est) <- paste("X",names(.est), sep=".")
  .ster = sqrt(diag(vcov(m))); names(.ster) <- paste("SE",names(.ster), sep=".")
  data.frame(c(.est, .ster)%>%t)
}

.compVarF <- function(df, err = 3e-2){
  err <- err*df$XSgl[1]
  ftr <- sum(df$FIDrvt^2)
  mutate(df, Wt = FIDrvt^2/ftr, WtT = Time*Wt, WtTT = Time^2*Wt)->df
  (sum(df$WtTT) - sum(df$WtT)^2)*ftr -> denom
  err/sqrt(denom)
}

.fit <- function(.sample, model, return.model = FALSE){
  w.g = rnorm(1, model@wFreq, .001)
  n.g = rnorm(1, model@Num0, .1*model@Num0)
  p.g = rnorm(1, model@Pol, .1*model@Pol)
  lam.decoh = model@decohLam
  phi = model@Phase
  
  cat("guessed frequency:", w.g, "\n\n")
  
  nls(Sgl ~ n*(1 + p*exp(lam.decoh*Time)*sin(frq*Time + phi)), data=.sample, start=list(n = n.g, p = p.g, frq = w.g)) -> m3
  
  if(return.model) return(m3) else return(.extract.stats(m3))
}

#### tests ####
#### SE(estimate(freq)) vs freq ####
varW0_test <- function(model, sampling, duration, wfreqs = c(.01,.1,.3,1,3,10)){
  require(doParallel); n.cores = detectCores()
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  rtns <- c(".extract.stats", ".fit", ".compVarF", "simSample","fiDer","expectation")
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
  
  stopCluster(clus)
  
  dat
  
}
