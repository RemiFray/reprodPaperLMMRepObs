
# sampler X ----

nimSamplerXmove <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 1
    S <- dim(model$latentObservation)[2]
    n <- control$n
    M <- control$M
    nM <- n+M
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    c <- nimSample(c(1:D))
    xProp <- model$x
    
    latentObsProp <- model$latentObservation
    latentIdxProp <- model$latentIndex
    
    log_qratio <- 0
    skip<-TRUE
    
    nu0 <- nimNumeric(S, value = 0.0)
    nu2 <- nimNumeric(S, value = 0.0)
    
    choice <- runif(1, 0, 1)
    
    ## ----------------------------------------------- add error ----
    if(choice <= 0.5){
      
      # now assuming its better to check if there is a 0 in the 
      #   one history sampled than sample one in which we know there is one
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
          # # calculate q(x|xprop) / x(xprop|x)
          X2 <- xProp[(n+1):nM] > 0
          qx <- 1/(sum(X2)*sum(nu2 == 2))
          qxprop <- 1/(sum(X0)*sum(nu0 == 0))
          log_qratio <- log(qx / qxprop)
          
          skip <- FALSE
        }
      }
    }
    
    # ----------------------------------------------- remove error ----
    else {
      
      X2 <- xProp[(n+1):nM] > 0
      if(sum(X2) != 0){
        # def nu2
        which_nu2 <- nimNumeric()
        which_nu2 <- which(X2)
        nu2i <- nimSample(which_nu2)+n
        nu2 <- latentObsProp[nu2i, ]
        # def t
        idxt <- nimNumeric() ; idxt <- which(nu2 == 2)
        t <- nimSample(idxt)
        # def nu1 and nu2
        nu1i <- which(latentIdxProp == 1+3^(t-1))[1]
        nu0 <- nu2 ; nu0[t] <- 0
        idxNu0i <- nimGetLatentIndex(nu0)
        if(any(latentIdxProp == idxNu0i))
          nu0i <- which(latentIdxProp == idxNu0i)[1]
        else{
          nu0i <- which(xProp[(n+1):nM] == 0)[1] + n # first free spot
          latentObsProp[nu0i, ] <- nu0
          latentIdxProp[nu0i] <- idxNu0i
        }
        # def xProp
        b <- numeric(nM)
        b[nu0i] <- + c
        b[nu1i] <- + c
        b[nu2i] <- - c
        xProp <- xProp + b
        if(xProp[nu2i] == 0)
          latentIdxProp[nu2i] <- 0
        # calculate q(x|xprop) / x(xprop|x)
        X0 <- xProp > 0
        qx <- 0.25/(sum(X0)*sum(nu0 == 0))
        qxprop <- 0.25/(sum(X2)*sum(nu2 == 2))
        log_qratio <- log(qx / qxprop)
        skip <- FALSE
      }
    }
    ## ------------------------------ metrolis ratio and acceptance ----
    if(all(xProp >= 0) & !skip){
      
      idxNu <- c(nu0i, nu1i, nu2i)
      model_lprop <- logRatio(idxNu, b, xProp, latentObsProp)
      
      log_MH_ratio <- model_lprop + log_qratio
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        model[[target]] <<- xProp
        model[["N"]] <<- sum(xProp)
        model[["latentObservation"]] <<- latentObsProp
        model[["latentIndex"]] <<- latentIdxProp
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N"), logProb = T)
      }
    }
  },
  methods = list(
    logRatio = function (idxNu=double(1), b=double(1), xProp = double(1), 
                         histories=double(2)) {
      N <- sum(model$x)
      Nprop <- N+sum(b)
      
      logRes <- lfactorial(Nprop)-lfactorial(N)
      for(i in 1:3){
        I <- idxNu[i]
        hist <- histories[I,]
        logRes <- logRes + lfactorial(model$x[I]) - lfactorial(xProp[I]) +
          b[I] * lProbHist(hist)
      }
      return(logRes)
      returnType(double(0))
    },
    lProbHist = function (hist=double(1)) {
      
      probHist <- 1 # pi_i
      for(t in 1:S){
        if(hist[t] == 0)
          probHist <- probHist * (1-model$capture[t])
        else if(hist[t] == 1)
          probHist <- probHist * model$capture[t] * model$alpha[1]
        else if(hist[t] == 2)
          probHist <- probHist * model$capture[t] * (1-model$alpha[1]) 
      }
      
      return(log(probHist))
      returnType(double(0))
    },
    reset = function () {})
)

# sampler unseen ----

nimSamplerX0 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    D <- if(!is.null(control$D))         control$D         else 10
    S <- dim(model$latentObservation)[2]
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    x0move <- nimSample(-D:D)
    xProp <- model$x
    if(xProp[1] + x0move >= 0){
      
      log_MH_ratio <- logRatio(x0move)
      
      u <- runif(1, 0, 1)
      if(u < exp(log_MH_ratio)){
        xProp[1] <- xProp[1] + x0move
        model[[target]] <<- xProp
        model[["D"]] <<- xProp[1]
        model[["N"]] <<- sum(xProp)
        
        copy(from = model, to = mvSaved, row = 1,
             nodes = c(calcNodes, "N", "D"), logProb = T)
      }
    }
  },
  methods = list(
    logRatio = function (x0move=double(0)) {
      N <- sum(model$x)
      Nprop <- N+x0move
      
      logRes <- lfactorial(Nprop)-lfactorial(N)
      
      hist <- nimNumeric(S)
      logRes <- logRes + lfactorial(model$x[1]) - lfactorial(model$x[1]+x0move) +
        x0move * lProbHist(hist)
      
      return(logRes)
      returnType(double(0))
    },
    lProbHist = function (hist=double(1)) {
      
      probHist <- 1 # pi_i
      for(t in 1:S){
        probHist <- probHist * (1-model$capture[t])
      }
      
      return(log(probHist))
      returnType(double(0))
    },
    reset = function () {} )
)

# sampler p ----

nimUpdatePt <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    S <- dim(model$latentObservation)[2]
    calcNodes <- model$getDependencies(target)
  },
  run = function() {
    
    apt <- nimNumeric(S, value = 1)
    bpt <- nimNumeric(S, value = 1)
    
    indexs <- which(model$x > 0)
    
    for (i in 1:length(indexs)){
      I <- indexs[i]
      for(t in 1:S){
        if(model$latentObservation[I,t] > 0)
          apt[t] <- apt[t] + model$x[I]# nb de captures total (> 0 dans les histo latents)
        else
          bpt[t] <- bpt[t] + model$x[I]# nb de non-capture total (0 dans les histo latents)
      }
    }
    
    new_p <- nimNumeric(S)
    for(t in 1:S){
      new_p[t] <- rbeta(1, apt[t],bpt[t])
    }
    
    
    model[[target]] <<- new_p
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( reset = function () {} )
)

# sampler alpha ----

nimUpdateA <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    prior <- control$prior
  },
  run = function() {
    
    aa = prior[1]
    ba = prior[2]
    
    indexs <- which(model$x > 0)
    
    for (i in 1:length(indexs)){
      I <- indexs[i]
      h <- model$latentObservation[I,]
      aa <- aa + model$x[I]*sum(h[] == 1)
      ba <- ba + model$x[I]*sum(h[] == 2)
    }
    new_a <- rbeta(1, aa, ba)
    
    model[[target]] <<- new_a
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = T)
  },
  methods = list( reset = function () {} )
)

