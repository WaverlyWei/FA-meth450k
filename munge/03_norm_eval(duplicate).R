################################################################################
# In this script we estimate the correlation between technical replicates as a
# way to evaluate the normalizations.
################################################################################

# Set up workspace
library(here)
library(limma)

# setting project and data directories
replicates_norm_data <- here("data", "replicates_norm_data.RData")

# only load data if replicates_norm_data does not exist
if (!file.exists(normalizePath(replicates_norm_data))) {

cor_est <- function(dat=pheno, exp){
# design setup
Condition <- factor(dat$dich_tce)
design <- model.matrix(~0+Condition)
colnames(design) <- levels(Condition)
# calculate correlation within subjects
corfit <- duplicateCorrelation(exp, design, block=dat$Subject)
return(corfit)
}

pheno_PS <- read.csv(file=here("data", "pheno_PS.csv"))
rownames(pheno_PS) <- pheno_PS[,1]
pheno_PS <- pheno_PS[,-1]
colnames(RGset_filtered) == rownames(pheno_PS)

load(here("data", "norm_data.RData"))

# Warning messages:
# 1: In glmgam.fit(dx, dy, coef.start = start, tol = tol, maxit = maxit,  :
#   Too much damping - convergence tolerance not achievable

# ENmix_BMIQ
corfit_ENmix_BMIQ <- cor_est(dat=pheno_PS, exp=mvals_ENmix_BMIQ)
corfit_ENmix_BMIQ$cor
# [1] 0.2720031
0.2933981

# ENmix_RCP
corfit_ENmix_RCP <- cor_est(dat=pheno_PS, exp=mvals_ENmix_RCP)
corfit_ENmix_RCP$cor
# [1] 0.2818666
0.3046645

# Funnorm
corfit_Funnorm <- cor_est(dat=pheno_PS, exp=mvals_Funnorm)
corfit_Funnorm$cor
# [1] 0.2768499
0.2696287

# Illumina
corfit_Illumina <- cor_est(dat=pheno_PS, exp=mvals_Illumina)
corfit_Illumina$cor
# [1] 0.3730618
0.4167187

# Quantile
corfit_Quantile <- cor_est(dat=pheno_PS, exp=mvals_Quantile)
corfit_Quantile$cor
# [1] 0.2554466
0.2692027

# Raw
corfit_Raw <- cor_est(dat=pheno_PS, exp=mvals_Raw)
corfit_Raw$cor
# [1] 0.241971
0.2610215

# Noob
corfit_Noob <- cor_est(dat=pheno_PS, exp=mvals_Noob)
corfit_Noob$cor
# [1] 0.2558961
0.2803224

# SWAN_ENmix
corfit_SWAN_ENmix <- cor_est(dat=pheno_PS, exp=mvals_SWAN_ENmix)
corfit_SWAN_ENmix$cor
# [1] 0.2630327
0.284417

# SWAN_Noob
corfit_SWAN_Noob <- cor_est(dat=pheno_PS, exp=mvals_SWAN_Noob)
corfit_SWAN_Noob$cor
# [1] 0.2441743
0.2628474

# transform back to the correlation scale and plot a histogram
# have the average correlation and the IQR

#############################################################################
###### save the relevant data
#############################################################################

save(
  corfit_ENmix_BMIQ, corfit_ENmix_RCP, corfit_Funnorm,
  corfit_Illumina,corfit_Quantile, corfit_Raw,
  corfit_Noob, corfit_SWAN_ENmix, corfit_SWAN_Noob,
  file = here("data", "replicates_norm_data.RData"))

} else {
 load(replicates_data)
}
