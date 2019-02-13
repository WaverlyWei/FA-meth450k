################################################################################
# Script for visualization, exploratory data analysis
################################################################################

library(here)
library(ChAMP)

pheno_PS <- read.csv(file=here("SuperFund/data", "pheno_PS.csv"))
rownames(pheno_PS) <- pheno_PS[,1]
pheno_PS <- pheno_PS[,-1]

run_svd <- function(betas, name, pheno){
  champ.SVD(beta = betas, rgSet=RGset_filtered, pd=pheno,
    RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("graphs", name))
}

load(here("data", "raw_data.RData"))
load(here("data", "norm_data.RData"))
load(here("data", "batchcorr_data.RData"))

pheno_sub <- pheno_PS[,-c(1,2,3,4,5,7,9,10,15,16,17,18,19,20,22,
  27,28,29,30,31,32,33,34,35)]
colnames(pheno_sub) <- c("Sample_Plate","TCE", "Sex", "Age","BMI","Smoking",
"Granulocyte","CD4", "CD8", "Bcell","NKcell")
pp <- pheno_sub[,c(8,9,10,11,7,5,4,3,6,2,1)]
run_svd(betas = betas_PS, pheno = pp, name = "SVD_1_afterFiltering")
run_svd(betas = betas_Illumina, pheno = pp,
  name = "SVD_2_afterNormalization")
betas_combat <- ENmix::M2B(mvals_Illumina_combat)
run_svd(betas = betas_combat, pheno = pp,
  name = "SVD_3_afterBatchCorr")
