###Color iPath map
rm(list=ls())
setwd("~/Box Sync/3rd Semester (May 2013 - November 2013)/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan")
require(RColorBrewer)

#Load
cntx.rxns   <- read.delim("rxns_cntx.txt",stringsAsFactors=F)
core.rxns   <- read.delim("rxns_core.txt",stringsAsFactors=F)
HMR2.rxnKEGG <- read.delim("../../data/cHMR3765_rxnKEGG.txt")

#Find KEGG id
cntx.iPath <- subset(HMR2.rxnKEGG,RXNID %in% cntx.rxns$RxnID,KEGG.ID)
core.iPath <- subset(HMR2.rxnKEGG,RXNID %in% core.rxns$RxnID,KEGG.ID)
cols <- tolower(brewer.pal(3,"RdBu")[c(1,3)])

#Prepare to write 2 columns id-color table
towrite <- data.frame(rbind(core.iPath,cntx.iPath),
                      c(rep(cols[1],nrow(core.iPath)),
                            rep(cols[2],nrow(cntx.iPath))))
towrite.f <- subset(towrite,KEGG.ID!="")
write.table(towrite.f,file = "corecntx2iPath.txt",sep = "\t",quote = F,row.names = F,col.names = F)
