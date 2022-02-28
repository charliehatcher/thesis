########################### script to diagnose MCMC convergence #################################

#load required packages 
library(coda)
library(data.table)
library(tidyverse)

setwd("/~/001_projects/002_selection/gctb_1.0_Linux/data/output/BayesNS/bw_MCMC")

list_files <- list.files("/~/001_projects/002_selection/gctb_1.0_Linux/data/output/BayesNS/bw_MCMC")

list <- list()

for (i in 1:length(list_files)){
  list[[i]] <- as.data.table(fread(list_files[i], header=T))
}

names(list) <- list_files
list2 <- lapply(list, raftery.diag)

############################ run on whole list of CpGs ###################################
results2 <- lapply(list2, function(x) {x <- (x[["resmatrix"]]); return(x)})
results3 <- as.data.frame(results2)
results4 <- t(results3)
results5 <- as.data.frame(results4[,c(1,5,8)])

results5$cpg <- rownames(results5)
results6<- results5 %>%
  filter(str_detect(cpg, '.I'))

#rename cpgs in cpg column 
results6$cpg <- gsub("*.mcmcsamples.Par.I", "", results6$cpg)
results6$cpg <- sub("^", "aa:",results6$cpg )

write.table(results6, "/~/MCMC_bw.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

rm(list=ls())

