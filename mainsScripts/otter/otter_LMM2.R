library(nimble)
library(tidyverse)
library(MCMCvis)
library(EnvStats) # normTrunc

# ------- Functions, samplers, distrib and model ----- ----

source("./functions/LMM2_Pt_DistribModel.R")
source("./functions/LMM2_Pt_functions.R")
source("./functions/LMM2_Pt_samplers.R")

# ------- Otter ------ ----

S <- 5
otter <- read.csv("data/otter.csv")


resDir <- paste0("./results/otter_LMM2")
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)



for(year in c("06", "07", "08", "10", "11", "12")){
  
  otterYear <- otter[grep(paste0("[0-9]+/[0-9]+/", year), otter$CollectionDate), ]
  
  dta_observed <- table(otterYear$OtterID, otterYear$CollectionDate) %>% 
    matrix(ncol = 5)

  # ------- Data ------ ----
  
  # MATRIX OF OBSERVED HISTORIES WITH SPACE FOR OTHER INDIV
  observation <- dta_observed[which(apply(dta_observed, 1, sum) > 0),]
  nb_obs <- nrow(observation)
  n <- 5*nb_obs
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
  for(j in 1:5){
    tmp <- addError(latObs2, errIndiv2, latIdx2, S)
    latObs2 <- tmp$latObs
    errIndiv2 <- tmp$errIndiv
    latIdx2 <- tmp$latIdx
  }
  
  # ------- Constructing nimble model ------ ----
  
  LMM2_Consts <- list(S = S, n = n)
  
  # Initialisation
  LMM2_Inits <- function(i) list(lambda = runif(S, 0.2, 0.8),
                                 alpha = runif(2, 0.6, 1)[i],
                                 N = c(sum(latentIndex>0), sum(latIdx2>0))[i],
                                 d = c(unseen_init, sum(latIdx2 == 1))[i],
                                 latentObservation = list(latentObservation, latObs2)[[i]],
                                 errIndiv = list(errIndiv, errIndiv2)[[i]],
                                 latentIndex = list(latentIndex, latIdx2)[[i]]
                                 )
  
  LMM2_Model <- nimbleModel(code = LMM2_Code, name = "LMM2", 
                              constants = LMM2_Consts,
                              data = list(), inits = LMM2_Inits(1))
  LMM2_Model$getLogProb("latentObservation")
  
  # Compilation modèle OK
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
  
  
  # Compilation MCMC
  LMM2_MCMC <- buildMCMC(LMM2_Conf)
  
  CLMM2_MCMC <- compileNimble(LMM2_MCMC, project = LMM2_Model,
                                showCompilerOutput = F, resetFunctions = T)
  
  
  # ------- run MCMC ------ ----
  
  inits <- list(LMM2_Inits(1), LMM2_Inits(2))
  
  burnin <- 1000
  nthin <- 10
  niter <- 10000+burnin

  samples <- runMCMC(CLMM2_MCMC, niter = niter, nburnin = burnin, 
                     thin = nthin, thin2 = 5,  nchains = 2,inits = inits, setSeed = c(777, 1234))
  
  
  save(samples, file = paste0(resDir, "/otter_LMM2_", year, ".Rdata"))
  
}


