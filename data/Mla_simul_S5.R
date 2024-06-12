library(nimble)
library(tidyverse)

# ----- Functions, samplers, distrib and model ----- ----

source("./functions/LMM2_Pt_functions.R")

# ----- Parameters ----- ----

S <- 5

l_test <- c(0.11, 0.23, 0.36, 0.51)
N_test <- c(500)
alpha_test <- c(0.8, 0.9, 0.95)

bias_plan <- data.frame(N = rep(N_test, each = length(l_test)*length(alpha_test)),
                        a = rep(rep(alpha_test, each = length(l_test)), length(N_test)),
                        l = rep(rep(l_test, length(alpha_test)), length(N_test)))

dtaDir <- paste0("./data/Mla_Simul_S",S)
if(!file.exists(dtaDir)) dir.create(dtaDir, recursive = TRUE)


# ------- Big loop ------ ----


for(line in 1:nrow(bias_plan)){
  set.seed(1234)
  
  lambda <- rep(bias_plan$l[line], S)
  alpha <- bias_plan$a[line]
  N <- bias_plan$N[line]
  
  for(iter in 1:10){
    
    # ------- Preparing data ------ ----
    
    capture_data <- simulCapture(N, S, lambda)
    res_observed <- simulIdentification(capture_data, alpha)
    
    save(res_observed,
         file = paste("./data/Mla_Simul_S",S, "/simul_S", S, "_N", N, "_a", alpha, 
                      "_l", lambda[1], "_iter", iter, ".Rdata", sep = ""))
  }
}

