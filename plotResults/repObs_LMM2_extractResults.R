# -------- Packages and functions -------- ----
library(dplyr)
library(MCMCvis)

wideSummary <- function(tmp){
  # if("D" %in% rownames(tmp)) idx <- c(2:4,6,8)
  # else idx <- c(1:3,5,7)
  
  rnames <- gsub("\\[|, |\\]", "", rownames(tmp), perl=TRUE)
  rnames <- gsub("capture", "p", rnames, perl=TRUE)
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


model <- "LMM2_Pt"

S <- c(5, 7, 9)
alpha <- c(0.8, 0.9, 0.95)
l_test <- c(0.11, 0.23, 0.36, 0.51)
N <- 500
prior <- c("noPrior")

# -------- Loop for extraction -------- ----

for(s in S){
  cat(s, "\n")
  testSummaries <- data.frame()
  for(l in l_test){
    cat("\t", l, "\n")
      for(a in alpha){
        cat("\t\t", a, "\n")
        params <- paste0("N", N, "_a", a, "_l", l)
        directory <- paste0("./results/simulations_S", s)
        files <- list.files(directory)
        for(iter in 1:10){
          
          fileExp <- paste0(params, "_iter", iter, ".Rdata")
          filePath <- files[which(grepl(fileExp, files))]
          load(file.path(directory, filePath))
          
          tmp <- data.frame(S=s, N=N, lambda=l, alpha=a, iter=iter) %>% 
            cbind(wideSummary(MCMCsummary(samples, 
                                          params = c("N", "alpha", "nbErrTot")))
            )
          testSummaries <- rbind(testSummaries, tmp)
        }
      }
    }

  write.csv(testSummaries,
            file=paste0("./results/repObs/", model, "_S", s, ".csv"), 
            row.names = F)
}


resS5 <- data.frame(model = rep("LMM2", 120)) %>% 
  cbind(read.csv(paste0("./results/repObs/", model, "_S5.csv")) )
resS7 <- data.frame(model = rep("LMM2", 120)) %>% 
  cbind(read.csv(paste0("./results/repObs/", model, "_S7.csv")) )
resS9 <- data.frame(model = rep("LMM2", 120)) %>% 
  cbind(read.csv(paste0("./results/repObs/", model, "_S9.csv")) )

allres <-  resS5 %>% 
  rbind(resS7) %>% 
  rbind(resS9)

write.csv(allres, file=paste0("./results/repObs/", model,
                              "_AllResults.csv"), row.names = F)
