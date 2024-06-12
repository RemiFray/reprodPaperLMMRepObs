
# ------- Distribution ------ ----

dLMM2 <- nimbleFunction(
  run = function(x = double(2), 
                 lambda = double(1),
                 alpha = double(0),
                 latentIndex = double(1),
                 errIndiv = double(2),
                 S = double(0), 
                 log = integer(0, default = 1) 
  ) {
    
    indexs <- which(latentIndex > 0)
    
    logProbData <- lfactorial(length(indexs))
    
    for(i in 1:length(indexs)){
      I <- indexs[i]
      hist <- x[I,]
      
      badid <- errIndiv[I, ]
      goodid <- hist - badid
      
      logProbHist <- 0
      for(t in 1:S){
        logProbHist <- logProbHist + log(dpois(hist[t], lambda[t]))
        if(hist[t] > 0)
          logProbHist <- logProbHist + log(dbinom(goodid[t], hist[t], alpha))
      }
      logProbData <- logProbData + logProbHist 
    }
    
    if(log) return(logProbData)
    return(exp(logProbData))
    returnType(double(0))
  })

rLMM2 <- nimbleFunction(
  run = function(n = integer(),
                 lambda = double(1),
                 alpha = double(0),
                 latentIndex = double(1),
                 errIndiv = double(2),
                 S = double(0)) {
    res <- matrix(nrow = length(latentIndex), ncol = S)
    returnType(double(2))
    return(res) ## dummy behavior
  })

deregisterDistributions("dLMM2")
registerDistributions(
  list(dLMM2 = list(BUGSdist = "dLMM2(lambda, alpha, latentIndex,
                                  errIndiv, S)",
                            types =c("value = double(2)", 
                                     "lambda = double(1)", 
                                     "latentIndex = double(1)", 
                                     "errIndiv = double(2)"
                                     ))))


# ------- Model Code ------ ----

LMM2_Code <- nimbleCode({
  
  # priors
  for(t in 1:S){
    lambda[t] ~ dgamma(shape = 1, rate = 1)
  }
  
  alpha ~ dbeta(1, 1)
  
  N ~ dunif(1, 10e4)
  D <- d
  
  for(t in 1:S){
    nbErr[t] <- sum(errIndiv[1:n, t])
  }
  nbErrTot <- sum(nbErr[1:S])
  
  # Likelihood, x
  latentObservation[1:n, 1:S] ~ dLMM2(lambda[1:S], alpha, latentIndex[1:n],
                                      errIndiv[1:n, 1:S], S)
  
})

