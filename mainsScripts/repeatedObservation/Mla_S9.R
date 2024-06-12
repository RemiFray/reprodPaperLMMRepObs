library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM2_Pt_DistribModel.R")
source("./functions/LMM2_Pt_functions.R")
source("./functions/LMM2_Pt_samplers.R")

# ----- Parameters ----- ----

modelType <- "LMMPoiss"
S <- 9

alpha_test <- c(0.8, 0.9, 0.95)
l_test <- c(0.11, 0.23, 0.36, 0.51)
N <- 500

priorType <- "uninformative"
priorAlpha <- c(1, 1)

bias_plan <- data.frame(a = rep(alpha_test, each = length(l_test)*10),
                        l = rep(rep(l_test, each = 10), length(alpha_test)),
                        iter = rep(rep(1:10, length(l_test)), length(alpha_test)))

burnin <- 100000
nthin <- 100
niter <- 10000*nthin+burnin


# ------- File management ------ ----

tmpDir <- paste0("tmp/", modelType, "_S", S, "_a", alpha, "_noPrior")
if(!file.exists(tmpDir)) dir.create(tmpDir, recursive = TRUE)

resDir <- paste0("./results/repObs/LMM2_S", S)
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)


# ------- Big loop ------ ----

for(line in 1:nrow(bias_plan)){
  
  # ------- Parameters ------ ----  
  N <- bias_plan$N[line]
  
  lambda <- rep(bias_plan$l[line], S)
  
  iter <- bias_plan$iter[line]
  
  # ------- Data ------ ----
  
  load(file = paste0("./data/Mla_Simul_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
                     "_l", lambda[1], "_iter", iter, ".Rdata"))
  
  dta_observed <- res_observed$observed
  
  # MATRIX OF OBSERVED HISTORIES WITH SPACE FOR OTHER INDIV
  observation <- dta_observed[which(apply(dta_observed, 1, sum) > 0),]
  nb_obs <- nrow(observation)
  n <- 5*N
  latentObservation <- matrix(nrow = n, ncol = S, data=0)
  latentObservation[1:nb_obs,] <- observation
  errIndiv <- matrix(nrow = n, ncol = S, data=0)
  latentIndex <- numeric(n)
  latentIndex[1:nb_obs] <- apply(observation, 1, nimGetLatentIndex, err = rep(0, S))
  unseen_init <- sample(5:50, 1)
  latentIndex[which(latentIndex == 0)[1:unseen_init]] <- 1
  rm(nb_obs)
  
  # DIFFERENT INITIALIZATION 
  set.seed(1248)
  latObs2 <- latentObservation
  errIndiv2 <- errIndiv
  latIdx2 <- latentIndex
  for(j in 1:40){
    tmp <- addError(latObs2, errIndiv2, latIdx2, S)
    latObs2 <- tmp$latObs
    errIndiv2 <- tmp$errIndiv
    latIdx2 <- tmp$latIdx
  }
  
  nb_err_init2 <- sum(errIndiv2>0)
  cat("Error added to second init: ", nb_err_init2, "\n")
  
  latIdx2[which(latIdx2 == 1)] <- 0
  latIdx2[which(latIdx2 == 0)[1:sample(5:sum(latIdx2>0), 1)]] <- 1
  
  # ------- Constructing nimble model ------ ----
  
  
  LMM2_Consts <- list(S = S, n = n)
  
  # Initialisation
  getlambdaInit <- function(i){
    latObs <- list(latentObservation, latObs2)[[i]]
    latIdx <- list(latentIndex, latIdx2)[[i]]
    latObs <- latObs[which(latIdx > 0), ]
    apply(latObs, 2, mean)
  }
  LMM2_Inits <- function(i){ list(lambda = getlambdaInit(i), 
                                 alpha = c(runif(2, 0.6, 1))[i],
                                 N = c(sum(latentIndex>0), sum(latIdx2>0))[i],
                                 d = c(unseen_init, sum(latIdx2 == 1))[i],
                                 latentObservation = list(latentObservation, 
                                                          latObs2)[[i]],
                                 errIndiv = list(errIndiv, errIndiv2)[[i]],
                                 latentIndex = list(latentIndex, latIdx2)[[i]]
  )}
  
  LMM2_Model <- nimbleModel(code = LMM2_Code, name = "LMM2", 
                            constants = LMM2_Consts,
                            data = list(), inits = LMM2_Inits(1))
  LMM2_Model$getLogProb("latentObservation")
  
  CLMM2_Model <- compileNimble(LMM2_Model, showCompilerOutput = FALSE)
  
  
  
  # ------- MCMC configuration ------ ----
  
  # Configuration of MCMC
  LMM2_Conf <- configureMCMC(LMM2_Model,
                             monitors = c("N", "D", "nbErrTot", "lambda", 
                                          "alpha"), 
                             thin = 1, thin2 = 5, 
                             print = TRUE)
  
  LMM2_Conf$removeSampler("lambda", "alpha", "N")
  LMM2_Conf$addSampler(target = "alpha", type = nimUpdateAlpha, 
                       control = list(prior = c(1, 1)))
  LMM2_Conf$addSampler(target = paste0("lambda[1:",S,"]"), type = nimUpdateLambda, 
                       control = list(priora = rep(1, S), priorb = 1))
  LMM2_Conf$addSampler(target = c("latentObservation"), type = nimSamplerXmove,
                       control = list(n = n, D = 1))
  LMM2_Conf$addSampler(target = c("latentObservation"), type = nimSamplerX0,
                       control = list(D = 5))
  LMM2_Conf$printSamplers()
  
  LMM2_MCMC <- buildMCMC(LMM2_Conf)
  
  CLMM2_MCMC <- compileNimble(LMM2_MCMC, project = LMM2_Model,
                              showCompilerOutput = F, resetFunctions = T)
  
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMM2_Inits(1), LMM2_Inits(2))
  
  samples <- runMCMC(CLMM2_MCMC, niter = niter, nburnin = burnin, 
                     thin = nthin, thin2 = 5,  nchains = 2,inits = inits, 
                     setSeed = c(777, 1234))

  save(samples,
       file = paste0(resDir,
                     "/LMM2_S", S, "_", priorType, "_N", N, "_a", alpha, 
                     "_l", lambda[1], "_iter", iter, ".Rdata"))
  
}



