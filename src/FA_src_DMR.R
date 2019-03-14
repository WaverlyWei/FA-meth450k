library(here)
source(here("src","utils.R"))
source(here("src","utils_annotation.R"))

pheno <- read.csv(here("data", "pheno_all.csv"))

sapply(pheno, class)
pheno <- run_class(pheno, fac = c(5,6,8,11:13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

design <- model.matrix(~fa+sex+age+infection+Gran_est+Mono_est+Bcell_est+
                         NK_est+CD4T_est+CD8T_est, data = pheno)

load(here("data", "combat.RData"))



#################################### DMR #######################################

## 1000 bootstraps
# run on terminal 
DMR <- run_DMR(mvals = mvals_combat, design = design, pheno = pheno)

save(DMR, file = here("data", "DMR.RData"),
     compress = TRUE)

DMR_sig <- DMR %>%
  filter(P.Value < .05)
dim(DMR_sig)

cat("DMR sig P")

DMR_sig <- DMR %>%
  filter(fwer < .05)
dim(DMR_sig)

cat("DMR sig fwer")