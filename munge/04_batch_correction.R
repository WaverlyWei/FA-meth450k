################################################################################
# In this script we use the normalized data sets and remove unwanted variation
# from batches (Sample_Plate) with ComBat.
################################################################################

# Set up workspace
library(sva)
library(here)

# setting data directory
batchcorr_data <- here("data", "batchcorr_data.RData")

# only load data from norm_data if batchcorr_data does not exist
if (!file.exists(normalizePath(batchcorr_data_data))) {

  print("batchcorr_data object not found; attempting to retrieve norm_data..")

  pheno_PS <- read.csv(file=here("SuperFund/data", "pheno_PS.csv"))
  rownames(pheno_PS) <- pheno_PS[,1]
  pheno_PS <- pheno_PS[,-1]

 ##############################################################################
 # ComBat
 # The ComBat function adjusts for known batches using an empirical Bayesian
 # framework. ComBat allows users to adjust for batch effects in datasets where
 # the batch cov is known, using methodology described in Johnson et al. 2007.
 ##############################################################################

 run_combat <- function(mvals, batch = pheno_PS$Sample_Plate) {

   mvals_combat <- ComBat(dat = mvals, batch)
   return(mvals_combat)
  }

# apply these functions to the normalizations

  load(here("data","norm_data.RData"))

  mvals_Illumina_combat <- run_combat(mvals=mvals_Illumina)

#############################################################################
###### save the relevant data
#############################################################################

save(mvals_Illumina_combat, file = here("SuperFund/data", "batchcorr_data.RData"))

} else {
 load(batchcorr_data)
}
