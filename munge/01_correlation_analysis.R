################################################################################
# In this  script we estimate the correlation between technical replicates
# as a way to choose a filtered dataset.
################################################################################

# Set up workspace
library(here)
library(limma)

# setting project and data directories
#replicates_data <- here("data", "replicates_data.RData")


# setting project and data directories
replicates_data <- here("SuperFund/data", "replicates_data.RData")


# only load data if replicates_data does not exist
if (!file.exists(normalizePath(replicates_data))) {
  
  ########################## load data ################################
  
  pheno_CGR <- read.csv(file=here("SuperFund/data", "pheno_CGR.csv"))
  pheno_P <- read.csv(file=here("SuperFund/data", "pheno_P.csv"))
  pheno_PS <- read.csv(file=here("SuperFund/data", "pheno_PS.csv"))
  
  load(here("SuperFund/data","raw_data.RData"))
  
  cor_est <- function(dat, exp){
    # design setup
    Condition <- factor(dat$fa)
    design <- model.matrix(~0+Condition)
    colnames(design) <- levels(Condition)
    # calculate correlation within subjects
    
    corfit <- duplicateCorrelation(exp, design, block=dat$Subject)
    return(corfit)
  }
  
  # CGR data
  corfit_CGR <- cor_est(dat=pheno_CGR, exp=ENmix::B2M(betas))
  corfit_CGR$cor
  ## ============== ERROR???? ============ #
  # [1] NaN
 
  # CGR data + additional probe filtering
  corfit_P <- cor_est(dat=pheno_P, exp=mvals_P)
  corfit_P$cor
  # [1] NaN
  
  # CGR data + additional probe and sample filtering
  corfit_PS <- cor_est(dat=pheno_PS, exp=mvals_PS)
  corfit_PS$cor
  # [1] NaN
  
  save(corfit_PS, corfit_P, corfit_CGR,
       file = here("data", "replicates_data.RData"))
  
} else {
  load(replicates_data)
}

  