################################################################################
# In this script we use the data that considers the selected probe & sample
# level filtering and perform several normalizations.
################################################################################

# Set up workspace
library(minfi)
library(ENmix)
library(here)

# setting data directory
norm_data <- here("data", "norm_data.RData")

# only load data from raw_data if norm_data does not exist
if (!file.exists(normalizePath(norm_data))) {

  print("norm_data object not found; attempting to retrieve raw_data..")

  load(here("SuperFund/data", "raw_data.RData"))
  pheno_PS <- read.csv(file=here("SuperFund/data", "pheno_PS.csv"))
  rownames(pheno_PS) <- pheno_PS[,1]
  pheno_PS <- pheno_PS[,-1]
  colnames(RGset_filtered) == rownames(pheno_PS)

  #############################################################################
  ###### normalization
  #############################################################################

  ################################# ENmix_BMIQ ################################

  Mset_processed <- preprocessENmix(RGset_filtered, exCpG=outCpG_PS)
  Mset_q1 <- norm.quantile(Mset_processed, method="quantile1")
  betas_ENmix_BMIQ <- bmiq.mc(Mset_q1)
  newBetas <- Harman::shiftBetas(betas_ENmix_BMIQ, shiftBy=1e-4)
  mvals_ENmix_BMIQ <- B2M(newBetas)
  betas_ENmix_BMIQ <- M2B(mvals_ENmix_BMIQ)
  dim(betas_ENmix_BMIQ)
  # [1] 399439    140

  ################################# ENmix_RCP #################################

  betas_ENmix_RCP <- rcp(Mset_q1)
  newBetas <- Harman::shiftBetas(betas_ENmix_RCP, shiftBy=1e-4)
  mvals_ENmix_RCP <- B2M(newBetas)
  betas_ENmix_RCP <- M2B(mvals_ENmix_RCP)
  dim(betas_ENmix_RCP)
  # [1] 399439    140

  ################################# Funnorm ###################################

  funnorm <- preprocessFunnorm(RGset_filtered, bgCorr = TRUE, dyeCorr = TRUE,
                   sex = ifelse(pheno_PS$sex == 1, "M",
                   ifelse(pheno_PS$sex == 0, "F", NA)), ratioConvert = FALSE)
  funnorm <- funnorm[!rownames(funnorm) %in% outCpG_PS]
  funnorm <- ratioConvert(funnorm, what = "both", keepCN = TRUE)
  funnorm <- mapToGenome(funnorm)
  betas <- assays(funnorm)$Beta
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_Funnorm <- B2M(newBetas)
  betas_Funnorm <- M2B(mvals_Funnorm)
  dim(betas_Funnorm)
  # [1] 399439    140

  ################################# Illumina ##################################

  illumina <- preprocessIllumina(RGset_filtered, bg.correct = TRUE,
    normalize = "controls")
  illumina <- illumina[!rownames(illumina) %in% outCpG_PS]
  illumina <- ratioConvert(illumina, what = "beta", keepCN = TRUE)
  illumina <- mapToGenome(illumina)
  ## =================== Question: What's mavls_combat??? =========== #
  assays(illumina)$M <- mvals_combat
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_Illumina <- B2M(newBetas)
  betas_Illumina <- M2B(mvals_Illumina)
  dim(betas_Illumina)
  # [1] 399439    140

  ################################# Quantile ##################################

  quan <- preprocessQuantile(RGset_filtered, fixOutliers = TRUE,
              quantileNormalize=TRUE, sex = ifelse(pheno_PS$sex == 1, "M",
              ifelse(pheno_PS$sex == 0, "F", NA)),
              stratified=TRUE)
  quan <- quan[!rownames(quan) %in% outCpG_PS]
  quan <- mapToGenome(quan)
  mvals <- assays(quan)$M
  betas <- M2B(mvals)
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  # No beta values were found to be 0 or 1. No shifts made.
  mvals_Quantile <- B2M(newBetas)
  betas_Quantile <- M2B(mvals_Quantile)
  dim(betas_Quantile)
  # [1]  403610     67

  ##################################### Raw ###################################

  raw <- preprocessRaw(RGset_filtered)
  raw <- raw[!rownames(raw) %in% outCpG_PS]
  raw <- ratioConvert(raw, what = "both", keepCN = TRUE)
  raw <- mapToGenome(raw)
  betas <- assays(raw)$Beta
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_Raw <- B2M(newBetas)
  betas_Raw <- M2B(mvals_Raw)
  dim(betas_Raw)
  # [1] 399439    140

  #################################### Noob ###################################

  Mset_noob <- preprocessNoob(RGset_filtered)
  noob <- Mset_noob[!rownames(Mset_noob) %in% outCpG_PS]
  noob <- ratioConvert(noob, what = "both", keepCN = TRUE)
  noob <- mapToGenome(noob)
  betas <- assays(noob)$Beta
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_Noob <- B2M(newBetas)
  betas_Noob <- M2B(mvals_Noob)
  dim(betas_Noob)
  # [1] 399439    140

  ################################ SWAN_ENmix #################################

  swan <- preprocessSWAN(RGset_filtered, mSet = Mset_processed, verbose = TRUE)
  swan <- swan[!rownames(swan) %in% outCpG_PS]
  swan <- ratioConvert(swan, what = "both", keepCN = TRUE)
  swan <- mapToGenome(swan)
  betas <- assays(swan)$Beta
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_SWAN_ENmix <- B2M(newBetas)
  betas_SWAN_ENmix <- M2B(mvals_SWAN_ENmix)
  dim(betas_SWAN_ENmix)
  # [1] 399439    140


  ################################ SWAN_Noob ##################################

  swan <- preprocessSWAN(RGset_filtered, mSet = Mset_noob, verbose = TRUE)
  swan <- swan[!rownames(swan) %in% outCpG_PS]
  swan <- ratioConvert(swan, what = "both", keepCN = TRUE)
  swan <- mapToGenome(swan)
  betas <- assays(swan)$Beta
  newBetas <- Harman::shiftBetas(betas, shiftBy=1e-4)
  mvals_SWAN_Noob <- B2M(newBetas)
  betas_SWAN_Noob <- M2B(mvals_SWAN_Noob)
  dim(betas_SWAN_Noob)
  # [1] 399439    140

  #############################################################################
  ###### save the relevant data
  #############################################################################

  ## Progress:
  # Figure out illumina part 
  # Done with quantile 
  # only save these two normalizations are fine 
  
  
  
  save(betas_ENmix_BMIQ, mvals_ENmix_BMIQ,
    betas_ENmix_RCP, mvals_ENmix_RCP,
    betas_Funnorm, mvals_Funnorm,
    betas_Illumina, mvals_Illumina,
    betas_Quantile, mvals_Quantile,
    betas_Raw, mvals_Raw,
    betas_Noob, mvals_Noob,
    betas_SWAN_ENmix, mvals_SWAN_ENmix,
    betas_SWAN_Noob, mvals_SWAN_Noob,
  file = here("data", "norm_data.RData"))

} else {
   load(norm_data)
}
