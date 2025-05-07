# -------- Packages and functions -------- ----
library(dplyr)
library(MCMCvis)

wideSummary <- function(tmp){
  rnames <- gsub("\\[|, |\\]", "", rownames(tmp), perl=TRUE)
  rnames <- gsub("capture", "p", rnames, perl=TRUE)
  rnames <- gsub("lambda", "l", rnames, perl=TRUE)
  rnames <- gsub("alpha", "a", rnames, perl=TRUE)
  cnames <- gsub("%", "", colnames(tmp), perl=TRUE)
  allCols <- as.character(sapply(rnames, function(x) paste(x, cnames, sep="_")))
  
  df <- as.matrix(tmp) %>% 
    t() %>% as.numeric() %>% 
    matrix(nrow=1) %>% as.data.frame()
  colnames(df) <- allCols
  df
}


# -------- Parameters -------- ----


folders <- c("otter_Yoshizaki", "otter_LMM2", "otter_Mt")
years <- c("06", "07", "08", "10", "11", "12")


for(i in 1:3){
  
  dta <- data.frame()
  files <- list.files(paste0("results/", folders[i]))
  model <- gsub("otter_", "", folders[i])
  
  for(j in 1:6){
    load(paste0("results/", folders[i], "/", files[j]))
    
    tmp <- data.frame(Model= model, 
                      year = years[j]) %>% 
      cbind(wideSummary(MCMCsummary(samples)))
    dta <- rbind(dta, tmp)
    
  }
  
  write.csv(dta,
            file=paste0("./results/otter_", model, "_Results.csv"), 
            row.names = F)
}


yoshi <- read.csv("./results/otter_Yoshizaki_Results.csv") %>% 
  select(Model, year, N = N_mean, N_2.5, N_97.5, N_sd, N_Rhat, N_n.eff,
         p1_mean, p1_2.5, p1_97.5, 
         p2_mean, p2_2.5, p2_97.5, 
         p3_mean, p3_2.5, p3_97.5, 
         p4_mean, p4_2.5, p4_97.5, 
         p5_mean, p5_2.5, p5_97.5) %>% 
  mutate(alpha = rep(NA, 6), alpha_sd = rep(NA, 6), nbErr = rep(NA, 6))
LMM2 <- read.csv("./results/otter_LMM2_Results.csv") %>% 
  select(Model, year, N = N_mean, N_2.5, N_97.5, N_sd, N_Rhat, N_n.eff,
         alpha = a_mean, alpha_sd = a_sd, 
         nbErr = nbErrTot_mean, 
         p1_mean = l1_mean, p1_2.5 = l1_2.5, p1_97.5 = l1_97.5, 
         p2_mean = l2_mean, p2_2.5 = l2_2.5, p2_97.5 = l2_97.5, 
         p3_mean = l3_mean, p3_2.5 = l3_2.5, p3_97.5 = l3_97.5, 
         p4_mean = l4_mean, p4_2.5 = l4_2.5, p4_97.5 = l4_97.5, 
         p5_mean = l5_mean, p5_2.5 = l5_2.5, p5_97.5 = l5_97.5)
Mt <- read.csv("./results/otter_Mt_Results.csv") %>% 
  select(Model, year, N = N_mean, N_2.5, N_97.5, N_sd, N_Rhat, N_n.eff,
         p1_mean, p1_2.5, p1_97.5, 
         p2_mean, p2_2.5, p2_97.5, 
         p3_mean, p3_2.5, p3_97.5, 
         p4_mean, p4_2.5, p4_97.5, 
         p5_mean, p5_2.5, p5_97.5) %>% 
  mutate(alpha = rep(NA, 6), alpha_sd = rep(NA, 6),  nbErr = rep(NA, 6))

AllModels <- rbind(yoshi, LMM2, Mt)
write.csv(AllModels,
          file=paste0("./results/otter_AllModels_Results.csv"), 
          row.names = F)
