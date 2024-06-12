library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/noSingle_DistribModel.R")

# ----- Parameters ----- ----

modelType <- "Yoshi"
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
nthin <- 10
niter <- 10000*nthin+burnin


# ------- File management ------ ----

resDir <- paste0("./results/repObs/", modelType, "_S",S)
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)


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
  dta_observed <- dta_observed[which(apply(dta_observed, 1, sum) > 1), ]
  dta_observed[which(dta_observed > 0)] <- 1
  dta_observed <- dta_observed[which(apply(dta_observed, 1, sum) > 1), ]
  
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
  M <- 0
  latentObservation <- matrix(nrow=n+M, ncol = S, data=0)
  latentObservation[1:n,] <- observation
  latentIndex <- numeric(n+M)
  latentIndex[1:n] <- summaryDta$index
  
  xInit <- numeric(n+M)
  xInit[1:n] <- summaryDta$Freq
  
  
  # ------- Constructing nimble model ------ ----
  
  
  noSingleConsts <- list(S = S, nbLatentObs=n+M,
                         latentObservation = latentObservation,
                         latentIndex = latentIndex)
  
  noSingleInits <- function() list(capture = rep(0.6, S))
  
  
  samples <- nimbleMCMC(code = noSingleCode, constants = noSingleConsts,
                         data = list(x = xInit), inits = noSingleInits(),
                         nchains = 2, niter = 10000, nburnin = 1000,
                         summary = FALSE, WAIC = FALSE, 
                         monitors = c('N','capture'))
  
  save(samples,
       file = paste0(resDir,
                     "/", modelType, "_S", S, "_", priorType, "_N", N, "_a", alpha, 
                     "_l", l, "_iter", iter, ".Rdata"))
  
}



