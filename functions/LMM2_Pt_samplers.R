
# x sampler ----
nimSamplerXmove <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    n <- length(model$latentIndex)
    D <- control$D
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    # ----------------------------------------------- Intro ----
    
    latentObsProp <- model$latentObservation
    latentIdxProp <- model$latentIndex
    errIndivProp <- model$errIndiv
    
    log_MH_ratio <- 0
    skip <- TRUE
    
    nu0 <- nimNumeric(S, value = 0.0)
    nu2 <- nimNumeric(S, value = 0.0)
    
    # ----
    for(nb_modif in 1:D){
    
    
    choice <- runif(1, 0, 1)
    # ----------------------------------------------- add error ----
    if(choice <= 0.5){
      # What time is available for adding an error
      available_t <- getAvailable0(latentIdxProp)
      if(length(available_t) > 0){
        skip <- FALSE
        # sample t ----
        t <- nimSample(available_t)
        # sample nu1 ----
        all_nu1i <- nimNumeric()
        all_nu1i <- which(latentIdxProp == 1+2^(t-1))
        nu1i <- all_nu1i[1]
        # sample nu0 ----
        all_nu0i <- nimNumeric()
        all_nu0i <- which(latentIdxProp > 0)
        nu0i <- nimSample(all_nu0i)
        while (nu0i == nu1i) {
          nu0i <- nimSample(all_nu0i)
        }
        nu0 <- latentObsProp[nu0i, ]
        nu2 <- nu0
        nu2[t] <- nu2[t] + 1
        
        # update nu0 as nu2 ----
        errIndivProp[nu0i, t] <- errIndivProp[nu0i, t] + 1 
        latentIdxProp[nu0i] <- nimGetLatentIndex(nu2, errIndivProp[nu0i, ])
        latentObsProp[nu0i, t] <- nu2[t]
        # put nu1 as unexisting ----
        latentIdxProp[nu1i] <- 0
        latentObsProp[nu1i, t] <- 0
        
        # q(x|xprop) / q(xprop|x) ----
        av2 <- getAvailable2(errIndivProp)
        all_nu2i <- nimNumeric()
        all_nu2i <- which(errIndivProp[, t] > 0)
        # [X | X'] / [X' | X]
        log_qx <- log(1/length(av2)) + log(1/length(all_nu2i))
        log_qxprop <- log(1/length(available_t)) + 
                      log(1/length(all_nu1i)) + log(1/(length(all_nu0i) - 1))
        log_qratio <- log_qx - log_qxprop
        # [Y | X', Z', N', ...] / [Y | X, Z, N, ...]
        N <- sum(model$latentIndex > 0)
        log_Lratio <- - log(N) + sum(model$lambda) - log(errIndivProp[nu0i, t]) +
          log(1 - model$alpha[1]) - log(model$alpha[1])
        # MH ratio
        log_MH_ratio <- log_MH_ratio + log_Lratio + log_qratio
      }
    }
    
    # ----------------------------------------------- remove error ----
    else {
      
      # What time is available for adding an error
      available_2 <- getAvailable2(errIndivProp)
      if(length(available_2) > 0){
        skip <- FALSE
        # sample t ----
        t <- nimSample(available_2)
        # sample nu2 ----
        all_nu2i <- nimNumeric()
        all_nu2i <- which(errIndivProp[, t] > 0)
        nu2i <- nimSample(all_nu2i)
        nu2 <- latentObsProp[nu2i, ]
        nu0 <- nu2
        nu0[t] <- nu0[t] - 1
        nu1i <- which(latentIdxProp == 0)[1]
        
        # put nu1 as seen ----
        latentIdxProp[nu1i] <- 1+2^(t-1)
        latentObsProp[nu1i, t] <- 1
        # update nu2 as nu0 ----
        errIndivProp[nu2i, t] <- errIndivProp[nu2i, t] - 1
        latentIdxProp[nu2i] <- nimGetLatentIndex(nu0, errIndivProp[nu2i, ])
        latentObsProp[nu2i, t] <- nu0[t]
        
        # q(x|xprop) / q(xprop|x) ----
        av0 <- getAvailable0(latentIdxProp)
        all_nu1i <- nimNumeric()
        all_nu1i <- which(latentIdxProp == 1+2^(t-1))
        all_nu0i <- nimNumeric()
        all_nu0i <- which(latentIdxProp > 0)
        possib0 <- length(all_nu0i - 1)
        # [X | X'] / [X' | X]
        log_qx <- log(1/length(av0)) + 
          log(1/length(all_nu1i)) + log(1/possib0)
        log_qxprop <- log(1/length(available_2)) + log(1/length(all_nu2i))
        log_qratio <- log_qx - log_qxprop
        # [Y | X', Z', N', ...] / [Y | X, Z, N, ...]
        Nprop <- sum(latentIdxProp > 0)
        log_Lratio <- log(Nprop)  - sum(model$lambda) + log(errIndivProp[nu2i, t]+1) +
          log(model$alpha[1]) - log(1 - model$alpha[1])
        # MH ratio
        log_MH_ratio <- log_MH_ratio + log_Lratio + log_qratio
      }
    }
    
    }
    # ------------------------------ metrolis ratio and acceptance ----
    
    if(!skip){
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        model[["N"]] <<- sum(latentIdxProp >  0)
        model[["D"]] <<- sum(latentIdxProp == 1)
        model[["latentObservation"]] <<- latentObsProp
        model[["latentIndex"]] <<- latentIdxProp
        model[["errIndiv"]] <<- errIndivProp
        
        model$calculate(c("nbErr", "nbErrTot"))
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N", "D", "errIndiv", "latentIndex", "nbErr"), logProb = T)
      }
    }
  },
  methods = list(
    getAvailable0 = function(latIdxProp = double(1)){
      available_t <- logical(S)
      for(t in 1:S){
        if(any(latIdxProp == 1+2^(t-1))) 
          available_t[t] <- TRUE
      }
      available_fully <- nimNumeric()
      available_fully <- which(available_t)
      return(available_fully)
      returnType(double(1))
    },
    getAvailable2 = function(errIndivProp = double(2)){
      available_t <- logical(S)
      for(t in 1:S){
        if(any( errIndivProp[, t] > 0)) 
          available_t[t] <- TRUE
      }
      available_fully <- nimNumeric()
      available_fully <- which(available_t)
      return(available_fully)
      returnType(double(1))
    },
    reset = function () {})
)



# Sampler unseen indiv ----

nimSamplerX0 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 10
    S <- dim(model$latentObservation)[2]
  },
  run = function() {
    
    x0move <- nimSample(c(-D:(-1), 1:D))
    if(model$D[1] + x0move >= 0){
      
      N <- sum(model$latentIndex > 0)
      Nprop <- N+x0move
      Dprop <- model$D[1] + x0move
      lprobHist <- - sum(model$lambda)
      # nb unseen binomial  N'! / n0'!(N'-n0')! pSeen^(N'-n0') pUnseen^n0'
      #                     N!  / n0!(N-n0)!    pSeen^(N-n0)   pUnseen^n0
      #               = (N'! / N) (n0! / n0'!) pUnseen^(n0'-n0)
      log_MH_ratio <- lfactorial(Nprop) - lfactorial(N) +
        lfactorial(model$D[1]) - lfactorial(Dprop) +
        x0move * lprobHist
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        latIdxProp <- model$latentIndex
        idx <- nimNumeric()
        nbMove <- abs(x0move)
        if(x0move > 0){
          idx <- which(latIdxProp == 0)[1:nbMove]
          latIdxProp[idx] <- 1
        }
        else{
          idx <- which(latIdxProp == 1)[1:nbMove]
          latIdxProp[idx] <- 0
        }
        
        model[["latentIndex"]] <<- latIdxProp
        model[["D"]] <<- Dprop
        model[["N"]] <<- Nprop
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c("latentIndex", "N", "D"), logProb = T)
      }
    }
  },
  methods = list(
    reset = function () {} )
)

# Sampler alpha ----
nimUpdateAlpha <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    prior <- control$prior
  },
  run = function() {
    
    aa <- prior[1]
    ba <- prior[2]
    
    nerr <- sum(model$errIndiv)
    ngood <- sum(model$latentObservation) - nerr
    
    aa <- aa + ngood
    ba <- ba + nerr
      
    new_a <- rbeta(1, aa, ba)
    
    model[[target]] <<- new_a
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( reset = function () {} )
)

# Sampler lambda ----
nimUpdateLambda <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    S <- length(model$lambda)
    priora <- control$priora
    priorb <- control$priorb
  },
  run = function() {
    
    a <- priora
    b <- priorb + model$N[1]
    
    new_lambda <- nimNumeric(S)
    for(t in 1:S){
      a[t] <- a[t] + sum(model$latentObservation[, t])
      new_lambda[t] <- rgamma(1, shape = a[t], rate = b)
    }
    
    model[[target]] <<- new_lambda
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( reset = function () {} )
)
