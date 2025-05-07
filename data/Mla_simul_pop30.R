library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM2_Pt_functions.R")

# ----- Parameters ----- ----

S <- 5

lambda <- rep(0.5, S)
N <- 30
alpha <- 0.95

# ------- Big loop ------ ----

set.seed(1234)

for(iter in 1:50){
  
  # ------- Preparing data ------ ----
  
  capture_data <- simulCapture(N, S, lambda)
  res_observed <- simulIdentification(capture_data, alpha)
  
  save(res_observed,
       file = paste("./data/Simul_pop30/simul_S", S, "_N", N, "_a", alpha, 
                    "_l", lambda[1], "_iter", iter, ".Rdata", sep = ""))
}
