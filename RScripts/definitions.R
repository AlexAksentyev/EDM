source("./RScripts/classes.R")
library(dplyr)

#### utility functions ####
.form <- function(x, n=3, spec="e") formatC(x,n,format=spec)
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

.compAnaWSE <- function(df, aerr = NA){
  if(is.na(aerr)) {cat("\t\t\t SE of the error is unavailable.\n\n"); return(NA)}
  ftr <- sum(df$FIDrvt^2)
  mutate(df, Wt = FIDrvt^2/ftr, WtT = Time*Wt)->df
  sum(df$WtT) -> MeanWT
  with(df, sum(Wt*(Time - MeanWT)^2))*ftr -> denom 
    # sum(Wt * (Time - MeanWT)^2) is analytically correct (I checked my derivation),
    # but the result is twice as large as the SE given by R.
    # sum(Wt* (Time - WtT)^2) doesn't seem to be correct analytically, but it seems to give R's SE
  aerr/sqrt(denom)
}

.fit <- function(.sample, model, return.model = FALSE){
  w.g = rnorm(1, model@wFreq, .01)
  n.g = rnorm(1, model@Num0, .1*model@Num0)
  p.g = rnorm(1, model@Pol, .1*model@Pol)
  lam.decoh = model@decohLam
  phi = model@Phase
  
  cat("Frequency guess:", w.g, "\n\n")
  
  nls(Sgl ~ n*(1 + p*exp(lam.decoh*Time)*sin(frq*Time + phi)), data=.sample, start=list(n = n.g, p = p.g, frq = w.g)) -> m3
  
  if(return.model) return(m3) else return(.extract.stats(m3)%>%cbind(SD.err = summary(m3)$sigma))
}

#### tests ####
#### SE(estimate(freq)) vs freq ####
varW0_test <- function(model, sampling, duration, wfreqs = c(.01,.1,.3,1,3,10)){
  require(doParallel); n.cores = detectCores()
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  rtns <- lsf.str(all=TRUE, envir=.GlobalEnv)
  clusterExport(clus, rtns)
  
  names(wfreqs) <- wfreqs
  llply(wfreqs, function(w) setValue(model, c("wFreq" = w))) %>%
    llply(
    function(.mod, .splstg, .Time){
      err = .splstg@rerror * .mod@Num0*.mod@Pol
      simSample(.splstg, .mod, .Time) -> .spl; .fit(.spl, .mod) %>%mutate(SEAN.frq=.compAnaWSE(.spl, err))-> .stats; 
      list("Stats" = .stats, "Sample" = .spl)
    }, .splstg = sampling, .Time = duration, 
    .parallel=TRUE, 
    .paropts = list(.export=c(rtns), .packages=c("dplyr"))
  ) -> dat
  
  stopCluster(clus)
  
  dat
  
}

#### total time test ####
.diagWrap <- function(model, samplings){
  require(doParallel); n.cores = detectCores()
  clus <- makeCluster(n.cores)
  registerDoParallel(clus)
  rtns <- lsf.str(all=TRUE, envir = .GlobalEnv) #c(".extract.stats", ".fit", ".compVarF", ".diagnosis", "simSample","fiDer","expectation")
  clusterExport(clus, rtns)
  
  tau = -1/model@decohLam; Ttot = 5*tau; sptt <- seq(5)*tau
  
  llply(samplings, function(s)
    simSample(s, model, Ttot) %>% mutate(Group = derivedFactor(
      "A" = Time <= sptt[1],
      "B" = Time <= sptt[2],
      "C" = Time <= sptt[3],
      "D" = Time <= sptt[4],
      "E" = Time <= sptt[5],
      .method = "first"
    )) %>% .diagnosis(model, s@rerror),
    .parallel = TRUE, 
    .paropts = list(.packages = c("dplyr", "mosaic"))
  ) -> dat; stopCluster(clus)
  
  dat
}
.diagnosis <- function(smpl, mod, rerr){
  N0P = mod@Num0*mod@Pol
  err = rerr * N0P
  smpl <- mutate(smpl, 
                 DDtNorm = mod@wFreq*FIDrvt/N0P,
                 Wt = FIDrvt^2/sum(FIDrvt^2)
  )
  
  smplslc <- smpl[seq(1,nrow(smpl), length.out=500),]
  psigscat <- ggplot(smplslc, aes(Time, Sgl, col=Group)) + geom_point() + 
    geom_line(
      aes(Time, XSgl), linetype=3, 
      data=data.frame(Group=NA, Time = seq(0, smpl[nrow(smpl), "Time"], length.out=500))%>%
        mutate(XSgl = expectation(mod, Time))
    ) + theme_bw()
  
  ddply(smpl, "Group", function(g) {
    mutate(g, Wt = FIDrvt^2/sum(FIDrvt^2)) ->g;
    .fit(g, mod)%>% mutate(SEAN.frq = .compAnaWSE(g, err), 
                 FItot = sum((g$FIDrvt/N0P)^2),
                 meanT = mean(g$Time),
                 WmeanT = sum(g$Time*g$Wt),
                 varT = var(g$Time),
                 WvarT = sum(g$Wt*(g$Time - WmeanT)^2)
          )}
  )  -> .stats
  .stats %>% transmute(Group, SE = SE.frq, SEAN = SEAN.frq, FItot) %>%
    reshape2::melt(id.vars=c("Group","FItot"), variable.name = "How", value.name="SE") -> .tstats
  pse <- ggplot(.tstats, aes(Group, SE, shape=How)) + geom_point() + theme_bw()
  
  pchar <- ggplot(.tstats) + geom_point(aes(FItot, SE, col=Group, shape=How)) +
    labs(x="Total information")
  
  pddt <- ggplot(smpl) + geom_density(aes(DDtNorm, col=Group)) +
    labs(x=expression(1/(N0P)~d/dt~s(w*t+phi)))
  
  list(
    "SglPlot" = psigscat, "SEPlot" = pse, 
    "CharPlot" = pchar, "ddtPlot" = pddt,
    "Sample" = smpl, "Stats" = .stats
  )
}

#### Compaction factor test ####
Comp_test <- function(model, samplings, Ttot){
  library(doParallel)
  makeCluster(detectCores()) -> clus; registerDoParallel(clus)
  rtns <- lsf.str(envir=.GlobalEnv, all=TRUE)
  clusterExport(clus, rtns)
  
  ldply(
    samplings,
    function(stm, model, Ttot){
      simSample(stm, model, Ttot) -> smpl
      .varWT(smpl)
    }, model, Ttot,
    .parallel = TRUE,
    .paropts = list(.packages = "dplyr"),
    .id = "Compact"
  ) -> dat; stopCluster(clus)
  
  dat
}

#### useful functions ####
## this is how I compute the denominator in .compAnaWSE
.varWT <- function(df){ 
  ftr = sum(df$FIDrvt^2)
  mutate(df, Wt = FIDrvt^2/ftr, WtT = Time*Wt)->df
  sum(df$WtT) -> MeanWT
  with(df, sum(Wt*(Time - MeanWT)^2)) -> VarWT
  
  data.frame("NUM" = nrow(df), "Ftr" = ftr, "MWT" = MeanWT, "VarT" = var(df$Time),"VarWT" = VarWT, Denom = ftr*VarWT)
}

## SingSam false formula
.SSSE <- function(sampling, model, Ttot){
  sderr = sampling@rerror*model@Num0*model@Pol
  w = model@wFreq; taud = -1/model@decohLam
  
  ftr = ifelse(model@decohLam == 0, pi/w/Ttot, (1-exp(-pi/w/taud))/(1-exp(-Ttot/taud)))
  sqrt(12 * ftr) *sderr/Ttot
}