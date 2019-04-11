library(here)
source(here("Desktop/SuperFund/src","utils.R"))

pheno <- read.csv(here("Desktop/SuperFund/data", "pheno_all.csv"))
sapply(pheno, class)
pheno <- run_class(pheno, fac = c(5,6,8,12,13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

# remove infection
design <- model.matrix(~fa+sex+age+Gran_est+Mono_est+Bcell_est+
                        NK_est+CD4T_est+CD8T_est, data = pheno)

# unadjusted design matrix
design_u <-  model.matrix(~fa, data = pheno)

load(here("Desktop/SuperFund/data", "combat.RData"))

# bonferroni threshold
0.05/nrow(mvals_combat)
# [1] 1.24277e-07


#################################### DMP #######################################
# Differential Mean DNA methylation 
DMP <- run_DMP(mvals = mvals_combat, design = design)
lambda(DMP[,5])
#[1] 1.154141


DMP_u <- run_DMP(mvals = mvals_combat, design = design_u)
lambda(DMP_u[,5])
#[1]  1.163215


save(DMP, file = here("Desktop/SuperFund/data", "DMP.RData"),
compress = TRUE)


save(DMP_u, file = here("Desktop/SuperFund/data", "DMP_unadjusted.RData"),
     compress = TRUE)


DMP_sig <- DMP %>%
filter(P.Value < .05)
dim(DMP_sig)
# [1] 26456     7


# unadjusted
DMP_sig_u <- DMP_u %>%
  filter(P.Value < .05)
dim(DMP_sig_u)
# [1] 24331     7


DMP_sig <- DMP %>%
filter(adj.P.Val < .05)
dim(DMP_sig)
# [1] 0 7

# unadjusted 
DMP_sig_u <- DMP_u %>%
  filter(adj.P.Val < .05)
dim(DMP_sig_u)
# [1] 0 7

#################################### DVP #######################################
# Differential Variable DNA methylation 
DVP <- run_DVP(mvals = mvals_combat, design = design)
lambda(DVP[,6])
# [1] 1.028389

# unadjusted
DVP_u <- run_DVP(mvals = mvals_combat, design = design_u)
lambda(DVP_u[,6])
# [1] 1.046189


save(DVP, file = here("Desktop/SuperFund/data", "DVP.RData"),
compress = TRUE)

# unadjusted
save(DVP_u, file = here("Desktop/SuperFund/data", "DVP_unadjusted.RData"),
     compress = TRUE)

DVP_sig <- DVP %>%
filter(P.Value < .05)
dim(DVP_sig)
# [1] 20408     7

# unadjusted 
DVP_sig_u <- DVP_u %>%
  filter(P.Value < .05)
dim(DVP_sig_u)
# [1] 21260     7



DVP_sig <- DVP %>%
filter(Adj.P.Value < .05)
dim(DVP_sig)
# [1] 10 7

# unadjusted
DVP_sig_u <- DVP_u %>%
  filter(Adj.P.Value < .05)
dim(DVP_sig_u)
# [1] 12 7

############ 
DVP_annot <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig)
write.csv(DVP_annot, here("Desktop/SuperFund/results","DVP.csv"),row.names = FALSE)

# unadjusted 
DVP_annot_u <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig_u)
write.csv(DVP_annot_u, here("Desktop/SuperFund/results","DVP_unadjusted.csv"),row.names = FALSE)



#################################### DMR #######################################
# DNA methylation across small regions
## 1000 bootstraps
# run on terminal 
DMR <- run_DMR(mvals = mvals_combat, design = design, pheno = pheno)

save(DMR, file = here("Desktop/SuperFund/data", "DMR.RData"),
compress = TRUE)

DMR_sig <- DMR %>%
filter(p.value < .05)
dim(DMR_sig)
# 183 14 

DMR_sig <- DMR %>%
filter(fwer < .05)
dim(DMR_sig)
#1 14 


write.csv(DMR_sig, here("Desktop/SuperFund/results","DMR.csv"),row.names = FALSE)

#################################### DMB #######################################
# large, cluster regions of CpG probes, blocks 
load(here("Desktop/SuperFund/data","normalized.RData"))
CN <- assays(normalized)$CN
filtered_CpG <- rownames(mvals_combat)
unfiltered_CpG <- rownames(CN)
outCpG <- unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]
CN <- CN[!rownames(CN) %in% outCpG,]
DMB <- run_DMB(mvals = mvals_combat, design = design, pheno = pheno,CN = CN)

save(DMB, file = here("Desktop/SuperFund/data", "DMB.RData"),
compress = TRUE)

DMB_sig <- DMB %>%
filter(p.value < .05)
dim(DMB_sig)
# 91 14

write.csv(DMB_sig, here("Desktop/SuperFund/results","DMB.csv"),row.names = FALSE)

DMB_sig <- DMB %>%
filter(fwer < .05)
dim(DMB_sig)
# 0 14 

load(file = here("Desktop/SuperFund/data", "DMP.RData"))
load(file = here("Desktop/SuperFund/data", "DVP.RData"))
load(file = here("Desktop/SuperFund/data", "DMB.RData"))
