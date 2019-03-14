library(tidyverse)

run_class <- function(pheno,fac = c(5,6,8,11:13), num = c(9,14:23)) {
  cols_all <- c(colnames(pheno))
  pheno[cols_all] <- lapply(pheno[cols_all], as.character)
  cols_num <- c(colnames(pheno[num]))
  pheno[cols_num] <- lapply(pheno[cols_num], as.numeric)
  cols_fac <- c(colnames(pheno[fac]))
  pheno[cols_fac] <- lapply(pheno[cols_fac], as.factor)
  return(pheno)
}
