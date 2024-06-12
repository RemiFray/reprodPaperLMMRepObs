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
  select(Model, year, N = N_mean, N_sd, p1 = p1_mean, p2 = p2_mean,
         p3 = p3_mean, p4 = p4_mean, p5 = p5_mean) %>% 
  mutate(alpha = rep(NA, 6), alpha_sd = rep(NA, 6), nbErr = rep(NA, 6))
LMM2 <- read.csv("./results/otter_LMM2_Results.csv") %>% 
  select(Model, year, N = N_mean, N_sd, alpha = a_mean, alpha_sd = a_sd, 
         nbErr = nbErrTot_mean, p1 = l1_mean, p2 = l2_mean,
         p3 = l3_mean, p4 = l4_mean, p5 = l5_mean)
Mt <- read.csv("./results/otter_Mt_Results.csv") %>% 
  select(Model, year, N = N_mean, N_sd, p1 = p1_mean, p2 = p2_mean,
         p3 = p3_mean, p4 = p4_mean, p5 = p5_mean) %>% 
  mutate(alpha = rep(NA, 6), alpha_sd = rep(NA, 6),  nbErr = rep(NA, 6))

AllModels <- rbind(yoshi, LMM2, Mt)
write.csv(AllModels,
          file=paste0("./results/otter_AllModels_Results.csv"), 
          row.names = F)
