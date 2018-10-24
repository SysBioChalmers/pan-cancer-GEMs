###Classify GEMs based on incl of core and cntx genes
rm(list=ls())
setwd("~/Documents/Academia/Projects/RxnLands/Codes/validateGEMs/")
library(ggplot2)
library(RColorBrewer)
source("performRFongenexpr.R")
source("performRFongenexpronlyT.R")

#1: Load data
##Gene expression and annotation
load("~/Documents/Data/annotation.RData")
genes <- unique(annotEntrez$ensembl_gene_id)

##Classifying genes: core only in T and cntx across T clusters
core.genes.onlyT  <- readLines("../analyzeGEMs/CompareToNormal/genes_core_onlyincancer.txt")
cntx.rxns.forcls  <- readLines("../analyzeGEMs/RxnsCoreVsPan/rxns_cntx_clust.txt")
###Find gene associated with cntx rxns
  HMR2.gpr   <- read.delim("../docs/cHMR3765_rxnGeneMatrix.txt")
  rownames(HMR2.gpr)<-HMR2.gpr[,1]
HMR2.gpr  <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]
cntx.genes.forcls <- c()
for (rxn in cntx.rxns.forcls){
  rxn.ex.genes     <- colnames(HMR2.gpr)[apply(HMR2.gpr[rxn,]!=0,2,any)]
  rxn.ex.genes     <- rxn.ex.genes[rxn.ex.genes %in% genes]
  cntx.genes.forcls <- c(cntx.genes.forcls,rxn.ex.genes)}
cntx.genes.forcls <- unique(cntx.genes.forcls)
###Translate to entrez
core.genes.onlyT  <- rownames(annotEntrez[annotEntrez$ensembl_gene_id%in%core.genes.onlyT,])
cntx.genes.forcls <- rownames(annotEntrez[annotEntrez$ensembl_gene_id%in%cntx.genes.forcls,])

#5: Regress label to rxn inclusion matrix
rf.cl <- performRFongenexpronlyT("Tclus",10)
rf.cl <- performRFongenexpr("Tclus",10)
rf.N  <- performRFongenexpr("TvsN",5)

#6: Compute multiclass AUC  
multiclass.roc(label.test,rf.pr)
