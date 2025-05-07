library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM_Pt_DistribModel.R")
source("./functions/LMM_Pt_functions.R")
source("./functions/LMM_Pt_samplers.R")

# ------- Otter ------ ----

S <- 5
otter <- read.csv("data/otter.csv")


resDir <- paste0("./results/otter_Mt")
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)



for(year in c("06", "07", "08", "10", "11", "12")){
  
  otterYear <- otter[grep(paste0("[0-9]+/[0-9]+/", year), otter$CollectionDate), ]
  
  dta_observed <- table(otterYear$OtterID, otterYear$CollectionDate) %>% 
    matrix(ncol = 5)
  dta_observed[which(dta_observed > 0)] <- 1
  summaryDta <- getSummaryCMR(dta_observed)
  
  if(!any(summaryDta$index == 1)) {
    summaryDta <- rbind(c(deparse(rep(0,S)), 1, 0), summaryDta)
    summaryDta$Freq <- as.numeric(summaryDta$Freq)
    summaryDta$index <- as.numeric(summaryDta$index)
  }
  
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
  
  # ------- Constructing nimble model ------ ----
  
  
  LMMPtConsts <- list(S = S, nbLatentObs=n+M)
  
  # Initialisation
  LMMPtInits <- function() list(capture = rep(runif(1, 0.1, 0.9), S),
                                alpha = 1, 
                                N = sum(xInit),
                                x = xInit,
                                latentObservation = latentObservation,
                                latentIndex = latentIndex)
  
  LMMPtModel <- nimbleModel(code = LMMPtCode, name = "Pt", 
                            constants = LMMPtConsts,
                            data = list(), inits = LMMPtInits())
  
  CLMMPt <- compileNimble(LMMPtModel, 
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
                              showCompilerOutput = F)
  
  # ------- run MCMC ------ ----
  set.seed(1234)
  inits <- list(LMMPtInits(), LMMPtInits(), LMMPtInits())
  
  burnin <- 1000
  nthin <- 1
  niter <- 10000+burnin
  
  system.time(
    samples <- runMCMC(CLMMPtMCMC, niter = niter, nburnin = burnin, 
                       thin = nthin, nchains = 3, inits = inits, setSeed = c(777, 1234, 555))
  )
  
  save(samples, file = paste0(resDir, "/otter_Mt_", year, ".Rdata"))
  
}



