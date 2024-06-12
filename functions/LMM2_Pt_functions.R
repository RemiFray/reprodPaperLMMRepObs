# ------- R functions ------ ----


# Simulation of data
simulCapture <- function(N, S, lambda){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    for(t in 1:S){
      cr[i,t] <- rpois(1, lambda) # Adding a succesful capture to cr
    }
  }
  
  return(cr)
  
}

# Simulation of data
simulCaptureNotPoiss <- function(N, S, n){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    for(t in 1:S){
      cr[i,t] <- sample(n, 1) # Adding a succesful capture to cr
    }
  }
  
  return(cr)
  
}

# Simulation of data
simulCaptureNegBin <- function(N, S, size, prob){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    for(t in 1:S){
      cr[i,t] <- rnbinom(1, size, prob) # Adding a succesful capture to cr
    }
  }
  
  return(cr)
  
}


# Simulation of data
simulIdentification <- function(captureData, alpha){
  
  N <- nrow(captureData)
  
  cr <- matrix(ncol = S, nrow = sum(captureData)+N, data = 0)
  cr[1:N, ] <- captureData
  errIndiv <- matrix(ncol = S, nrow = N, data = 0)
  
  n_err_made <- 0
  
  for(i in 1:N){
    for(t in 1:S){
      if(cr[i, t] > 0){
        for(it in 1:cr[i, t]){
          if(runif(1) > alpha){
            n_err_made <- n_err_made + 1
            n <- N + n_err_made
            cr[i, t] <- cr[i, t] - 1
            cr[n, t] <- 1
            errIndiv[i, t] <- errIndiv[i, t] + 1 
          }
        }
      }
    }
    
  }
  end <- N + n_err_made
  cr <- cr[1:end, ]
  
  return(list(observed = cr, errIndiv = errIndiv))
  
}

# Adds error to generate different xInit
addError <- function(latentObservation, errIndiv, latentIndex, S){
  
  skip <- TRUE
  errIndivProp <- errIndiv
  latentObsProp <- latentObservation
  latIdxProp <- latentIndex
  
  available_t <- logical(S)
  for(t in 1:S){
    if(any(latIdxProp == 1+2^(t-1))) available_t[t] <- TRUE
  }
  available_t <- which(available_t)
  if(length(available_t) > 0){
    skip <- FALSE
    t <- nimSample(available_t)
    
    all_nu1i <- nimNumeric()
    all_nu1i <- which(latIdxProp == 1+2^(t-1))
    nu1i <- all_nu1i[1]
    
    all_nu0i <- nimNumeric()
    all_nu0i <- which(latIdxProp > 0)
    nu0i <- nimSample(all_nu0i)
    while (nu0i == nu1i) {
      nu0i <- nimSample(all_nu0i)
    }
    nu0 <- latentObsProp[nu0i, ]
    nu2 <- nu0
    nu2[t] <- nu2[t] + 1
    
    errIndivProp[nu0i, t] <- errIndivProp[nu0i, t] + 1 
    latIdxProp[nu0i] <- nimGetLatentIndex(nu2, errIndivProp[nu0i, ])
    latentObsProp[nu0i, t] <- nu2[t]
    
    latIdxProp[nu1i] <- 0
    latentObsProp[nu1i, t] <- 0
  }
  
  if(!skip) 
    return(list(latObs=latentObsProp, latIdx=latIdxProp,
                errIndiv=errIndivProp, err=1))
  else return(list(latObs=latentObservation, latIdx=latentIndex,
                   errIndiv=errIndiv, err=0))
}


# ------ Functions ------ ------


# sample uniform double function compiled
nimSample <-  nimbleFunction(
  run = function(x = double(1)){
    l <- length(x)
    ind <- rcat(1, nimRep(1/l, l))
    nsamp <- x[ind]
    
    returnType(double(0))
    return(nsamp)
  })
# CnimSample <- compileNimble(nimSample)

nimGetLatentIndex <- nimbleFunction(
  run = function(h = double(1), err = double(1)){
    if(sum(h) == 0) res <- 1
    else if(sum(h) == 1 & sum(err) == 0) {
      t <- which(h == 1)
      res <- 1 + 2^(t[1] - 1)
    }
    else res <- 0.5
    return(res)
    returnType(double(0))
  })
# CnimGetLatentIndex <- compileNimble(nimGetLatentIndex)

