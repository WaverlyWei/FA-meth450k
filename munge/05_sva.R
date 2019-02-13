################################################################################
# In this script we use the normalized data sets and remove unwanted variation
# with surrogate variable analysis (SVA).
################################################################################

# Set up workspace
library(sva)
library(here)

# setting data directory
design_sva <- here("SuperFund/data", "design_sva.RData")

# only load data from norm_data if cellcorr_data does not exist
if (!file.exists(normalizePath(design_sva))) {

  print("design_sva object not found; retrieving batchcorr_data..")

  pheno_PS <- read.csv(file=here("SuperFund/data", "pheno_PS.csv"))
  rownames(pheno_PS) <- pheno_PS[,1]
  pheno_PS <- pheno_PS[,-1]

  # drop <- dplyr::filter(pheno_PS,TCE_M_C == 1)
  # dropList <- drop$Target_ID
  # pheno <- pheno_PS[!rownames(pheno_PS) %in% dropList,]

  ##############################################################################
  # SVA
  # The goal of the sva is to remove all unwanted sources of variation while
  # protecting the contrasts due to the primary variables included in mod.
  # This leads to the identification of features that are consistently different
  # between groups, removing all common sources of latent variation.
  ##############################################################################

  run_sva <- function(mvals, var) {

    # mod is the the model matrix being used to fit the data

    # The full model includes terms for both the adjustment variables and
    # the variables of interest. The variable of interest might an indicator of
    # cancer versus control. The adjustment variables could be the age of the
    # patients, the sex of the patients, and a variable like the date the
    # arrays were processed.

    mod <- model.matrix(~var, data=pheno_PS)

#model.matrix(~dich_tce+Lymphocyte+Monocyte+Granulocyte
#system is computationally singular: reciprocal condition number = 4.27583e-17
#model.matrix(~dich_tce+Lymphocyte+Granulocyte
#system is computationally singular: reciprocal condition number = 4.96905e-17

    # mod0 is the null model being compared when fitting the data

    # The null model is a model matrix that includes terms for all of the
    # adjustment variables but not the variables of interest.

    mod0 <- model.matrix(~1, data=pheno_PS)

    sva <- sva(mvals, mod, mod0,  n.sv = 9)

# https://www.biostars.org/p/121489/
# get a "cleaned" version of the matrix with the surrogate variables removed.
    # cleanY <- function(y, mod, svs) {
    #   X <- cbind(mod, svs)
    #   Hat <- solve(t(X) %*% X) %*% t(X)
    #   beta <- (Hat %*% t(y))
    #   rm(Hat)
    #   gc()
    #   P <- ncol(mod)
    #   return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
    # }

   design_sva <- cbind(mod, sva$sv)
   return(design_sva)
 }

  load(here("data","batchcorr_data.RData"))

  mvals_combat <- mvals_Illumina_combat
  # mvals_combat <- mvals_combat[,!colnames(mvals_combat) %in% dropList]
  mod0 <- model.matrix(~1, data=pheno_PS)
  #design_sva_fa <- run_sva(mvals = mvals_combat, var = fa)
  
  mod_fa <- model.matrix(~fa, data=pheno_PS)
  sva_fa <- sva(mvals_combat, mod_fa, mod0,  n.sv = 9)
  design_sva_fa <- cbind(mod_fa, sva_fa$sv)
  
  #design_sva_C <- run_sva(mvals = mvals_combat, var = TCE_M_C)
  
  mod_ppm <- model.matrix(~log(fa_ppm), data=pheno_PS)
  sva_ppm <- sva(mvals_combat, mod_ppm, mod0,  n.sv = 9)
  #design_sva_ppm <- run_sva(mvals = mvals_combat, var = log(fa_ppm))
  design_sva_ppm <- cbind(mod_ppm, sva_ppm$sv)
  
  #save(design_SV, file = here("data", "design_smk.RData"))

  #############################################################################
  ###### save the relevant data
  #############################################################################

  save(design_sva_fa, design_sva_ppm,
    file = here("SuperFund/data", "design_sva.RData"))

  } else {
   load(design_sva)
  }
