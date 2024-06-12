

dLMMPt <- nimbleFunction(
  run = function(x = double(1), 
                 capture = double(1),
                 identification = double(0),
                 latentObservation = double(2),
                 latentIndex = double(1),
                 S = double(0), 
                 log = integer(0, default = 1) 
  ) {
    
    logProbData <- lfactorial(sum(x)) # full likelihood
    indexs <- which(x > 0)
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- latentObservation[I,]
      logProbData <- logProbData - lfactorial(x[I]) # log( 1/Prod(x_i!) )
      probHist <- 1 # pi_i
      for(t in 1:S){
        if(hist[t] == 0)
          probHist <- probHist * (1-capture[t])
        else if(hist[t] == 1)
          probHist <- probHist * capture[t] * identification
        else if(hist[t] == 2)
          probHist <- probHist * capture[t] * (1-identification) 
      }
      logProbData <- logProbData + x[I] * log(probHist) # x_i * log( pi_i ))
    }
    
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  })

rLMMPt <- nimbleFunction(
  run = function(n = integer(), capture = double(1), 
                 identification = double(0), latentObservation = double(2),
                 latentIndex = double(1), S = double(0)) {
    returnType(double(1))
    return(latentIndex) ## dummy behavior
  })

deregisterDistributions("dLMMPt")
registerDistributions(
  list(dLMMPt = list(BUGSdist = "dLMMPt(capture, identification, 
                latentObservation, latentIndex, S)",
                            types =c("value = double(1)", "capture = double(1)", 
                                     "latentObservation = double(2)", 
                                     "latentIndex = double(1)"))))


# ------- Model Code ------ ----

LMMPtCode <- nimbleCode({
  
  # priors
  for(t in 1:S){
    capture[t] ~ dbeta(1.0, 1.0)
  }
  
  alpha ~ dbeta(1.0, 1.0)
  
  # N ~ dnegbin(size=30,p=30/530)
  N ~ dunif(1, 10e4)
  D <- x[1]
  
  # Likelihood, x
  x[1:nbLatentObs] ~ dLMMPt(capture[1:S], alpha, 
                            latentObservation[1:nbLatentObs, 1:S], 
                            latentIndex[1:nbLatentObs], S)
  
})

