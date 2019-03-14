##############################################################################
# ComBat
# The ComBat function adjusts for known batches using an empirical Bayesian
# framework. ComBat allows users to adjust for batch effects in datasets where
# the batch cov is known, using methodology described in Johnson et al. 2007.
##############################################################################
library(sva)
library(here)

source(here("SuperFund/munge","utils.R"))

pheno <- read.csv(here("SuperFund/data", "pheno_all.csv"))
sapply(pheno, class)
pheno <- run_class(pheno, fac = c(5,6,8,11:13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

load(here("SuperFund/data", "filtered.RData"))

batch <- pheno$Sample_Plate
mod <- model.matrix(~fa+sex+age+infection+Gran_est+Mono_est+Bcell_est+
                      NK_est+CD4T_est+CD8T_est, data=pheno)


mvals_combat <- ComBat(dat = mvals, batch, mod)

save(mvals_combat,file = here("SuperFund/data", "combat.RData"),
     compress = TRUE)
