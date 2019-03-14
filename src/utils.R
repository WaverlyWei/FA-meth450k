library(limma)
library(minfi)
library(missMethyl)
library(doParallel)

run_DMP <- function(mvals, design){
  lFit <- limma::lmFit(object = mvals, design = design)
  eFit <- limma::eBayes(lFit, robust = TRUE)
  pTop <- limma::topTable(eFit, coef = 2, num = Inf, sort.by = "p")
  pTop <- pTop[order(pTop$P.Value), , drop = FALSE]
  ProbeID <- rownames(pTop)
  DMP <- data.frame(ProbeID, pTop)
  return(DMP)
}

run_DVP <- function(mvals, design){
  fitvar <- varFit(mvals, design = design, coef = c(1,2))
  topDV <- topVar(fitvar, coef=2, number = 399439)
  topDV <- data.frame(Probe_ID = rownames(topDV), topDV)
  return(topDV)
}

run_DMR <- function(mvals, design, pheno){
  registerDoParallel(cores = 1)
  GRset <- makeGenomicRatioSetFromMatrix(mat = mvals,
    pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M")
  Bumphunter <- minfi::bumphunter(object = GRset, design = design, coef = 2,
    B = 1000, type = "M", pickCutoff=TRUE, nullMethod="bootstrap",
    maxGap=1000, smooth = TRUE, smoothFunction = loessByCluster)
  dmr_results <- Bumphunter$table
  return(dmr_results)
}

run_DMB <- function(mvals, design, pheno, CN){
  registerDoParallel(cores = 1)
  GRset <- makeGenomicRatioSetFromMatrix(mat = mvals,
    pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M")
  assays(GRset)$CN <- CN
  cluster <- cpgCollapse(GRset, what = "M")
  blocks <- minfi::blockFinder(cluster$object, design, coef = 2,
  what = "M", nullMethod = "bootstrap", smooth = TRUE, cluster = NULL,
  pickCutoff=TRUE, B = 1000, smoothFunction = loessByCluster)
  dmb_results <- blocks$table
  return(dmb_results)
}

lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE),
na.rm=TRUE) / qchisq(0.5, df=1)

## ============== Where to get annotation data ============ ##

annotate_probes <- function(probeID_column, sig_dat){
  load(here("SuperFund/data", "450k_annot.RData"))
  sig_probes <- as.character(sig_dat[,probeID_column])
  sig_dat$Probe_ID <- sig_dat[,probeID_column]
  annot <- annotation[(annotation$Probe_ID %in% sig_probes),]
  dat <- merge(sig_dat, annot, by = "Probe_ID")
  return(dat)
}



######## run_class ########

library(tidyverse)
run_class <- function(pheno,fac = c(5,6,8,11:13), num = c(9,14:23)) {
  cols_all <- c(colnames(pheno))
  pheno[cols_all] <- lapply(pheno[cols_all], as.character)
  cols_num <- c(colnames(pheno[num]))
  pheno[cols_num] <- lapply(pheno[cols_num], as.numeric)
  cols_fac <- c(colnames(pheno[fac]))
  pheno[cols_fac] <- lapply(pheno[cols_fac], as.factor)
  return(pheno)
}
