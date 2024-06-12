library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM_Pt_DistribModel.R")
source("./functions/LMM_Pt_functions.R")
source("./functions/LMM_Pt_samplers.R")

# ----- Parameters ----- ----

modelType <- "Mt"
S <- 7

alpha_test <- c(0.8, 0.9, 0.95)
l_test <- c(0.11, 0.23, 0.36, 0.51)
N <- 500

priorType <- "uninformative"
priorAlpha <- c(1, 1)

bias_plan <- data.frame(a = rep(alpha_test, each = length(l_test)*10),
                        l = rep(rep(l_test, each = 10), length(alpha_test)),
                        iter = rep(rep(1:10, length(l_test)), length(alpha_test)))

burnin <- 10000
nthin <- 100
niter <- 10000*nthin+burnin


# ------- File management ------ ----

tmpDir <- paste0("tmp/", modelType, "_S", S, "_a", alpha, "_noPrior")
if(!file.exists(tmpDir)) dir.create(tmpDir)

resDir <- paste0("./results/repObs/", modelType, "_S",S)
if(!file.exists(resDir)) dir.create(resDir)


# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  
  # ------- Parameters ------ ----  
  
  alpha <- bias_plan$a[line]
  l <- bias_plan$l[line]
  
  iter <- bias_plan$iter[line]
  
  # ------- Data ------ ----
  
  load(file = paste0("./data/Mla_Simul_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
                     "_l", l, "_iter", iter, ".Rdata"))
  
  dta_observed <- res_observed$observed
  dta_observed[which(dta_observed > 0)] <- 1
  summaryDta <- getSummaryCMR(dta_observed)
  
  # MATRIX OF ORDERED OBSERVED HISTORIES
  observation <- as.data.frame(t(sapply(summaryDta$history, 
                                        function(h) eval(parse(text=h)))))
  rownames(observation) <- NULL
  observation$index <- apply(observation, 1, nimGetLatentIndex)
  observation <- observation %>% 
    arrange(index) %>% select(-index) %>% as.matrix()
  colnames(observation) <- NULL
  
  # MATRIX OF HISTORIES AND SPACE FOR OTHERS
  n <- nrow(observation)
  M <- min(sum(summaryDta$Freq[-1]), 3^S-n)
  latentObservation <- matrix(nrow=n+M, ncol = S, data=0)
  latentObservation[1:n,] <- observation
  latentIndex <- numeric(n+M)
  latentIndex[1:n] <- summaryDta$index
  
  xInit <- numeric(n+M)
  xInit[1:n] <- summaryDta$Freq
  
  # MATRIX A (just to check if xInit if valid)
  A <- getA(latentObservation, observation)
  
  # INITIAL X TO START MCMC
  set.seed(1248)
  xInit2 <- xInit
  latObs2 <- latentObservation
  latIdx2 <- latentIndex
  for(j in 1:40){
    tmp <- addError(xInit2, latObs2, latIdx2, S, n, M, n+M)
    xInit2 <- tmp$x
    latObs2 <- tmp$latObs
    latIdx2 <- tmp$latIdx
    
    A2 <- getA(latObs2, observation)
    if(!(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2)))
      print(j)
  }
  cat(all(summaryDta$Freq[-1] == A %*% xInit), "  ")
  A2 <- getA(latObs2, observation)
  cat(all(xInit2 >= 0) & all(summaryDta$Freq[-1] == A2 %*% xInit2), "  ")
  print(sum(xInit-xInit2))
  
  # ------- Constructing nimble model ------ ----
  
  
  LMMPtConsts <- list(S = S, nbLatentObs=n+M)
  
  # Initialisation
  LMMPtInits <- function(i) list(capture = rep(0.5, S),
                                 alpha = 1, 
                                 N = c(sum(xInit),sum(xInit2))[[i]],
                                 x = list(xInit, xInit2)[[i]],
                                 latentObservation = list(latentObservation, latObs2)[[i]],
                                 latentIndex = list(latentIndex, latIdx2)[[i]])
  
  LMMPtModel <- nimbleModel(code = LMMPtCode, name = "Pt", 
                            constants = LMMPtConsts,
                            data = list(), inits = LMMPtInits(1))
  
  CLMMPt <- compileNimble(LMMPtModel, 
                          dirName = tmpDir,
                          showCompilerOutput = FALSE)
  
  
  # ------- MCMC configuration ------ ----
  
  # Configuration of MCMC
  LMMPtConf <- configureMCMC(LMMPtModel,
                             monitors = c("N", "D", "capture"), 
                             thin = 1,
                             print = TRUE)
  
  LMMPtConf$removeSampler("x", "capture", "alpha", "N")
  LMMPtConf$addSampler(target = paste0("capture[1:",S,"]"),
                       type = nimUpdatePt)
  LMMPtConf$addSampler(target = c("x"), type = nimSamplerX0,
                       control = list(D = 5))
  LMMPtConf$printSamplers()
  
  LMMPtMCMC <- buildMCMC(LMMPtConf)
  
  
  CLMMPtMCMC <- compileNimble(LMMPtMCMC, project = LMMPtModel,
                              dirName = tmpDir,
                              showCompilerOutput = F)
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMMPtInits(1), LMMPtInits(1))
  
  system.time(
    samples <- runMCMC(CLMMPtMCMC, niter = niter, nburnin = burnin, 
                       thin = nthin, nchains = 2,inits = inits, setSeed = c(777, 1234))
  )
  
  save(samples,
       file = paste0(resDir,
                     "/", modelType, "_S", S, "_", priorType, "_N", N, "_a", alpha, 
                     "_l", l, "_iter", iter, ".Rdata"))
  
}



