library(here)
source(here("SuperFund/src","utils.R"))

pheno <- read.csv(here("SuperFund/data", "pheno_all.csv"))
sapply(pheno, class)
pheno <- run_class(pheno, fac = c(5,6,8,11:13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

design <- model.matrix(~fa+sex+age+infection+Gran_est+Mono_est+Bcell_est+
                        NK_est+CD4T_est+CD8T_est, data = pheno)

load(here("data", "combat_Shanghai.RData"))

# bonferroni threshold
0.05/nrow(mvals_combat)
# [1] 1.24277e-07
#################################### DMP #######################################
DMP <- run_DMP(mvals = mvals_combat, design = design)
lambda(DMP[,5])
#[1] 1.154677

save(DMP, file = here("SuperFund/data", "DMP.RData"),
compress = TRUE)

DMP_sig <- DMP %>%
filter(P.Value < .05)
dim(DMP_sig)
# [1] 26516     7

DMP_sig <- DMP %>%
filter(adj.P.Val < .05)
dim(DMP_sig)
# [1] 0 7

#################################### DVP #######################################
DVP <- run_DVP(mvals = mvals_combat, design = design)
lambda(DVP[,6])
# [1] 1.027794

save(DVP, file = here("SuperFund/data", "DVP.RData"),
compress = TRUE)

DVP_sig <- DVP %>%
filter(P.Value < .05)
dim(DVP_sig)
# [1] 20393     7

DVP_sig <- DVP %>%
filter(Adj.P.Value < .05)
dim(DVP_sig)
# [1] 10 7

############ 
DVP_annot <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig)
write.csv(DVP_annot, here("SuperFund/results","DVP.csv"),row.names = FALSE)



#################################### DMR #######################################

## 1000 bootstraps
# run on terminal 
DMR <- run_DMR(mvals = mvals_combat, design = design, pheno = pheno)

save(DMR, file = here("SuperFund/data", "DMR_Shanghai.RData"),
compress = TRUE)

DMR_sig <- DMR %>%
filter(P.Value < .05)
dim(DMR_sig)
#

DMR_sig <- DMR %>%
filter(fwer < .05)
dim(DMR_sig)
#

#################################### DMB #######################################
DMB <- run_DMB(mvals = mvals_s, design = design_s, pheno = pheno_s)

save(DMB, file = here("data", "DMB.RData"),
compress = TRUE)

DMB_sig <- DMB %>%
filter(P.Value < .05)
dim(DMB_sig)
#

DMB_sig <- DMB %>%
filter(fwer < .05)
dim(DMB_sig)
#
