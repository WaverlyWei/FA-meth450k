library(here)
library(ChAMP)
library(minfi)
library(ENmix)

source(here("SuperFund/munge", "utils.R"))

pheno <- read.csv(here("SuperFund/data", "pheno_all.csv"))
sapply(pheno, class)
pheno <- run_class(pheno,  fac = c(5,6,8,11:13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

## 
load(here("SuperFund/data", "raw_data.RData"))
load(here("SuperFund/data", "normalized.RData"))
# crucial that the methylation data columns and phenotype rows correspond

# Check matching
colnames(RGset) == rownames(pheno)
colnames(assays(normalized)$M) == rownames(pheno)
colnames(assays(normalized)$Beta) == rownames(pheno)

# detP - so we can filter the CpGs that didn't successfully hybridize to array
detP <- detectionP(RGset)
det <- detP[match(rownames(assays(normalized)$Beta), rownames(detP)),]
rownames(assays(normalized)$Beta) == rownames(det)

filtered <- champ.filter(beta=assays(normalized)$Beta,
                         M=assays(normalized)$M, pd = pheno, detP = det)

# Filtering probes with a detection p-value above 0.01.
# Removing 14187 probes.
# If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
# 
# Filtering NoCG Start
# Only Keep CpGs, removing 2661 probes from the analysis.
# 
# Filtering SNPs Start
# Using general 450K SNP list for filtering.
# Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
# Removing 56907 probes from the analysis.
# 
# Filtering MultiHit Start
# Filtering probes that align to multiple locations as identified in Nordlund et al
# Removing 11 probes from the analysis.
# 
# Filtering XY Start
# Filtering probes located on X,Y chromosome, removing 9419 probes from the analysis.
# 
# Updating PD file
# 
# Fixing Outliers Start
# Replacing all value smaller/equal to 0 with smallest positive value.
# Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]
# 
# All filterings are Done, now you have 402327 probes and 72 samples.


betas <- filtered$beta
mvals <- filtered$M

save(betas, mvals,file=here("SuperFund/data","filtered.RData"),
     compress = TRUE)

champ.QC(beta = betas, pheno= pheno$fa, mdsPlot=TRUE,
         densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,
         Feature.sel="None", resultsDir=here::here("SuperFund/graphs"))
champ.SVD(beta = betas, rgSet=RGset, pd=pheno, RGEffect=TRUE,
          Rplot=FALSE, resultsDir=here::here("SuperFund/graphs", "SVD_"))
