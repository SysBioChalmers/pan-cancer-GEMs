###Classify GEMs based on incl of core and cntx rxns
rm(list=ls())
setwd("~/Documents/Academia/Projects/RxnLands/Codes/validateGEMs/")

source("performRFonrxninclmat.R")

#Load data
##Gene expression
load("../data/preprocessed.RData") 
##Model inclusion matrices for T
  rxnIncMatFile     <- c("../analyzeGEMs/rxnInclMat.txt")
  rxnInclMat        <- read.delim(rxnIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat  <- rxnInclMat[,-grep("X.1",colnames(rxnInclMat))]
##Model inclusion matrices for T
  rxnIncMatFile_N     <- c("rxnIncMatr_N.txt")
  rxnInclMat_N        <- read.delim(rxnIncMatFile_N,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat_N        <- rxnInclMat_N[,-grep("X.1",colnames(rxnInclMat_N))]
##Rxn annot in HMR2
HMR2.rxnTable <- read.delim("../data/cHMR3765_rxnTable.txt",stringsAsFactors=F)
  rownames(HMR2.rxnTable) <- HMR2.rxnTable$RxnID
##Classifying rxns: core only in T and cntx across T clusters
core.rxns.onlyT   <- readLines("../analyzeGEMs/CompareToNormal/rxns_core_onlyincancer.txt")
cntx.rxns.forcls  <- readLines("../analyzeGEMs/RxnsCoreVsPan/rxns_cntx_clust.txt")

#2: Create a factor
  PSID         <- rownames(y$samples)
  PSID.s       <- gsub(pattern="-",replacement=".",x=PSID)
  filter       <- PSID.s%in%colnames(rxnInclMat)
  mapper       <- match(colnames(rxnInclMat),PSID.s[filter])
fGlobalCT.m  <- fGlobalCT[filter][mapper]
  names(fGlobalCT.m)  <- PSID.s[filter][mapper]

#5: Regress label to rxn inclusion matrix
res  <- performRFonrxninclmat("Tclus",4)
print(res[["rf"]])
print(res[["multi.roc"]])
