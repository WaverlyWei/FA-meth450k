# prep data for plotting probes

library(ENmix)
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

# load mvals_combat
betas <- M2B(mvals_combat)

data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
oth <- data.frame(Other)
loc <- data.frame(Locations)
oth$Probe_ID <- rownames(Other)
dusp22 <- dplyr::filter(oth, grepl('DUSP22', UCSC_RefGene_Name))
rownames(dusp22) <- dusp22$Probe_ID
betas_dusp22 <- betas[rownames(betas) %in% dusp22$Probe_ID,]
tt <- merge(betas_dusp22, dusp22, by = "row.names")
rownames(tt) <- tt$Row.names
loc_dusp22 <- loc[rownames(loc) %in% rownames(tt),]
t1 <- merge(tt, loc_dusp22, by = "row.names")
rownames(t1) <- rownames(tt)
write.csv(t1, here("Desktop/SuperFund/results", "dusp22.csv"))

pheno <- read.csv(here("Desktop/SuperFund/data", "pheno_all.csv"))
ph <- t(pheno)
colnames(ph) <- ph[3,]

# stop here
# exposure 
pp <- data.frame(ph[6,])
p1 <- data.frame(t(pp))

# skip subsetting 
#t2 <- t1[,c(3:142)]
#rownames(t2) <- rownames(t1)
#t3 <- base::rbind(t2,p1)

#  
t1 <- t1[,-1]
names <- colnames(p1)
# merged_dat add FA as a row, merged by TargetID
merged_dat <- merge(t1,p1,by = names, all = TRUE)
rownames(merged_dat) <- c(t1$Row.names, "FA")
# order dat by "pos", increasing order 
ordered_dat <- merged_dat[order(merged_dat$pos),]

write.csv(merged_dat, here("Desktop/SuperFund/results", "dusp22_pheno.csv"))
write.csv(ordered_dat, here("Desktop/SuperFund/results", "ordered_dusp22_pheno.csv"))

# read in 
dusp22 <- read.csv(here("Desktop/SuperFund/results", "ordered_dusp22_pheno.csv"))

