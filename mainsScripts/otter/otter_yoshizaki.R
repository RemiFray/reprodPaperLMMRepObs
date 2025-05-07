library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/noSingle_DistribModel.R")

# ------- Otter ------ ----

S <- 5
otter <- read.csv("data/otter.csv")


resDir <- paste0("./results/otter_Yoshizaki")
if(!file.exists(resDir)) dir.create(resDir, recursive = TRUE)



for(year in c("06", "07", "08", "10", "11", "12")){
  
  otterYear <- otter[grep(paste0("[0-9]+/[0-9]+/", year), otter$CollectionDate), ]
  
  dta_observed <- table(otterYear$OtterID, otterYear$CollectionDate) %>% 
    matrix(ncol = 5)
  dta_observed[which(dta_observed > 0)] <- 1
  dta_observed <- dta_observed[which(!(apply(dta_observed, 1, sum) == 1)), ]
  
  summaryDta <- getSummaryCMR(dta_observed)
  
  
  # ------- Preparation ------ ----
  
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
  
  noSingleInits <- function() list(capture = rep(runif(1, 0.5, 0.9), S))
  
  
  samples <- nimbleMCMC(code = noSingleCode, constants = noSingleConsts,
                        data = list(x = xInit), inits = noSingleInits(),
                        nchains = 3, niter = 10000, nburnin = 1000,
                        summary = FALSE, WAIC = FALSE, 
                        monitors = c('N','capture'))
  
  save(samples, file = paste0(resDir, "/otter_", year, ".Rdata"))
  
  
}
