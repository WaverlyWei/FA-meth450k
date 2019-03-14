################################################################################
# normalize the data and estimate the counts for each cell type
################################################################################

library(here)
library(minfi)
library(ENmix)
library(tidyverse)
source(here("SuperFund/munge","utils.R"))

load(here("data", "raw_data.RData"))
pheno <- read.csv(here("data", "pheno_PS.csv"))

# Add cell counts columns
cell_counts <- estimateCellCounts(RGset,
  referencePlatform = "IlluminaHumanMethylation450k", returnAll = TRUE)

colnames(cell_counts$counts) <- c("CD8T_est","CD4T_est","NK_est","Bcell_est",
                           "Mono_est","Gran_est")
# create blood cell dataframe
bc <- data.frame(cell_counts$counts)
bc$Target_ID <- rownames(bc)
ph <- left_join(pheno, bc)
# check variable class
sapply(ph,class)

# factors: sample_plate, sample_well, sample_group, study, sex
# alc, inf, smoke, fa, grp_comb

# numeric: age, all blood cell counts 

ph <- run_class(ph,  fac = c(5,6,8,11:13), num = c(9,14:23))
write.csv(ph, here("SuperFund/data", "pheno_all.csv"), row.names = FALSE)

normalized <- cell_counts$normalizedData
bb <- M2B(assays(normalized)$M)
assays(normalized)$Beta <- bb
save(normalized, file = here("SuperFund/data", "normalized.RData"),
     compress = TRUE)
