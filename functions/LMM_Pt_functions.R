# ------- R functions ------ ----

latentObsGenerates <- function(histo, vectHists){
  index2 <- which(histo == 2)
  hist2 <- histo
  hist2[index2] <- 0
  
  a <- integer(length(vectHists))
  
  if(any(hist2 > 0)){
    a[which(vectHists == deparse(hist2))] <- 1
  }
  
  if(length(index2) > 0){
    for(i in index2){
      histErr <- numeric(length(histo))
      histErr[i] <- 1
      a[which(vectHists == deparse(histErr))] <- 1
    }
  }
  
  a
}

getA <- function(latentObservation, observation){
  vectHists <- apply(observation, 1, deparse)
  A <- t(apply(latentObservation, 1, latentObsGenerates, 
               vectHists = vectHists))
  t(A)[-which(vectHists==deparse(numeric(S))),]
}


# Simulation of data
simulCMR <- function(N, S, pcapture, identification){
  
  cr <- matrix(ncol = S, nrow = N, data = 0) # Matrix of capture NbIndiv * NbSample
  
  for(i in 1:N){
    for(t in 1:S){
      if(runif(1) < pcapture[t]) { #capture
        if (runif(1) < identification)
          cr[i,t] <- 1 # Adding a succesful capture to cr
        else{
          ccr <- cr
          cr <- matrix(nrow = nrow(ccr)+1, ncol = S, data = 0)
          cr[1:nrow(ccr),] <- ccr[,]
          cr[nrow(cr),t] <- 1 # add one hist with unique observation
        }
      }
    }
  }
  
  return(cr)
  
}


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

# get the x without errors that verify y = A*x
getXinit <- function(latentObservation, dtaCR, init){
  xInit <- apply(latentObservation, 1, function(x) sum(x != 2))
  xInit[xInit != S] <- 0
  xInit[xInit > 0] <- dtaCR$Freq
  xInit[1] <- init
  xInit
}


# Adds error to generate different xInit
addError <- function(x, latentObservation, latentIndex, S, n, M, nM){
  
  c <- 1
  xProp <- x
  
  latentObsProp <- latentObservation
  latentIdxProp <- latentIndex
  
  nu0 <- nimNumeric(S, value = 0.0)
  nu2 <- nimNumeric(S, value = 0.0)
  
  X0 <- xProp > 0
  # def nu0
  which_nu0 <- nimNumeric()
  which_nu0 <- which(X0)
  nu0i <- nimSample(which_nu0)
  nu0 <- latentObsProp[nu0i, ]
  idxt <- nimNumeric() ; idxt <- which(nu0 == 0)
  if(length(idxt) > 0){
    t <- nimSample(idxt)
    # def nu1 and nu2
    nu1iIdx <- which(latentIdxProp == 1+3^(t-1))
    if(length(nu1iIdx) > 0){
      nu1i <- nu1iIdx[1]
      nu2 <- nu0 ; nu2[t] <- 2
      idxNu2 <- nimGetLatentIndex(nu2)
      # nu2 already exists
      if(any(latentIdxProp==idxNu2))
        nu2i <- which(latentIdxProp==idxNu2)[1]
      # nu2 needs to be added
      else {
        nu2i <- which(xProp[(n+1):nM] == 0)[1] + n # first free spot
        latentObsProp[nu2i, ] <- nu2
        latentIdxProp[nu2i] <- idxNu2
      }
      b <- numeric(nM)
      b[nu0i] <- - c
      b[nu1i] <- - c
      b[nu2i] <- + c
      xProp <- xProp + b
    }
  }
  
  
  if(all(xProp >= 0)) 
    return(list(x=xProp, latObs=latentObsProp, latIdx=latentIdxProp))
  else return(list(x=x, latObs=latentObservation, latIdx=latentIndex))
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
  run = function(h = double(1)){
    s <- 1 
    for(t in 1:length(h)) s <- s + h[t]*3^(t-1)
    return(s)
    returnType(double(0))
  })
# CnimGetLatentIndex <- compileNimble(nimGetLatentIndex)



