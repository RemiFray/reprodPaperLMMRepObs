# -------- Packages and functions -------- ----
library(dplyr)
library(MCMCvis)

wideSummary <- function(tmp){
  
  rnames <- gsub("\\[|, |\\]", "", rownames(tmp), perl=TRUE)
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


model <- "repObs_pop30"

S <- 5
alpha <- 0.95
l_test <- 0.5
N <- 30
prior <- c("noPrior")

# -------- Loop for extraction -------- ----

for(s in S){
  testSummaries <- data.frame()
  for(l in l_test){
    for(a in alpha){
      params <- paste0("N", N, "_a", a, "_l", l)
      directory <- "./results/LMM2_pop30"
      files <- list.files(directory)
      for(iter in 1:50){
        
        cat(iter, "  ")
        fileExp <- paste0(params, "_iter", iter, ".Rdata")
        filePath <- files[which(grepl(fileExp, files))]
        load(file.path(directory, filePath))
        
        tmp <- data.frame(S=s, N=N, lambda=l, alpha=a, iter=iter) %>% 
          cbind(wideSummary(MCMCsummary(samples, 
                                        params = c("N", "alpha", "lambda", "nbErrTot")))
          )
        testSummaries <- rbind(testSummaries, tmp)
      }
    }
  }
  
  write.csv(testSummaries,
            file=paste0("./results/", model, ".csv"), 
            row.names = F)
}

