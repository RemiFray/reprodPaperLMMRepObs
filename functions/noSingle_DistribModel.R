

dnoSingle <- nimbleFunction(
  run = function(x = double(1), 
                 capture = double(1),
                 tauStar = double(0),
                 latentObservation = double(2),
                 latentIndex = double(1),
                 S = double(0), 
                 log = integer(0, default = 1) 
  ) {
    
    logProbData <- lfactorial(sum(x))
    indexs <- which(x > 0)
    
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- latentObservation[I,]
      logProbData <- logProbData - lfactorial(x[I])
      probHist <- 1 
      for(t in 1:S){
        if(hist[t] == 0)
          probHist <- probHist * (1-capture[t])
        else if(hist[t] == 1)
          probHist <- probHist * capture[t]
      }
      logProbData <- logProbData + x[I] * log(probHist / tauStar)
    }
    
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  })

rnoSingle <- nimbleFunction(
  run = function(n = integer(), capture = double(1), 
                 tauStar = double(0), latentObservation = double(2),
                 latentIndex = double(1), S = double(0)) {
    returnType(double(1))
    return(latentIndex) ## dummy behavior
  })

deregisterDistributions("dnoSingle")
registerDistributions(
  list(dnoSingle = list(BUGSdist = "dnoSingle(capture, tauStar, latentObservation, 
                                    latentIndex, S)",
                            types =c("value = double(1)", "capture = double(1)", 
                                     "latentObservation = double(2)", 
                                     "latentIndex = double(1)"))))


# ------- Model Code ------ ----

noSingleCode <- nimbleCode({
  
  # priors
  for(t in 1:S){
    capture[t] ~ dbeta(1.0, 1.0)
  }
  
  tau[S+1] <- prod(1 - capture[1:S])
  for(t in 1:S) tau[t] <- capture[t] * prod(1 - capture[1:S]) / (1-capture[t])
  
  tauStar <- 1 - sum(tau[1:(S+1)])
  
  N <- sum(x[1:nbLatentObs]) / tauStar
  
  # Likelihood, x
  x[1:nbLatentObs] ~ dnoSingle(capture[1:S],  tauStar,
                            latentObservation[1:nbLatentObs, 1:S], 
                            latentIndex[1:nbLatentObs], S)
  
})



# ------ Functions ------ ------

# Counts of histories in sorted vector
getSummaryCMR <- function(dta){
  summaryDta <- as.data.frame(table(apply(dta, 1, deparse)), 
                              stringsAsFactors = F)
  summaryDta$index <- sapply(summaryDta$Var1, 
                             function(h) nimGetLatentIndex(eval(parse(text=h))))
  summaryDta <- summaryDta %>% 
    arrange(index) %>% 
    select(history=Var1, index, Freq)
  summaryDta
}

nimGetLatentIndex <- nimbleFunction(
  run = function(h = double(1)){
    s <- 1 
    for(t in 1:length(h)) s <- s + h[t]*3^(t-1)
    return(s)
    returnType(double(0))
  })