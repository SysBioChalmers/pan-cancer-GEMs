###Define core rxns
rm(list=ls())
library(gplots)
library(ggplot2)
setwd("~/Box Documents/3rd Semester (May 2013 - November 2013)/Projects/MutLands/ProjectCodes/MutLands/analyzeGEMs/ModelComparison/RxnsCoreVsPan/")
percSampleWRxns <- 0.99

#Load files
load("../../../data/preprocessed_AllMutCID.RData") 
rxnIncMatFile <- c("../../rxnInclMat.txt")#,"../geneInclMat.txt","../metInclMat.txt")
rxnInclMat    <- read.delim(rxnIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat    <- rxnInclMat[,-grep("X.1",colnames(rxnInclMat))]
geneIncMatFile <- c("../../geneInclMat.txt")#,"../geneInclMat.txt","../metInclMat.txt")
geneInclMat    <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
geneInclMat    <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
HMR3676.rxnTable <- read.delim("../../../docs/cHMR3765_rxnTable.txt")
HMR3675.gpr  <- read.delim("../../../docs/cHMR3765_rxnGeneMatrix.txt")
rownames(HMR3675.gpr)<-HMR3675.gpr[,1]
HMR3675.gpr <- HMR3675.gpr[,-c(1,ncol(HMR3675.gpr))]

# Create a factor
globalPSID.s <- gsub(pattern="-",replacement=".",x=globalPSID)
filter       <- globalPSID.s%in%colnames(rxnInclMat)
mapper       <- match(colnames(rxnInclMat),globalPSID.s[filter])
fGlobalCT.f  <- fGlobalCT[filter]
fGlobalCT.m  <- fGlobalCT.f[mapper]

#Find n° rxns per model
GEMS.rxnXPSID          <- apply(rxnInclMat==1,2,sum)
GEMS.rxnXPSID.sort.ind <- sort.int(GEMS.rxnXPSID,index.return = T)
GEMS.rxnXPSID.sort     <- GEMS.rxnXPSID[GEMS.rxnXPSID.sort.ind$ix]
barplot(GEMS.rxnXPSID.sort,names.arg = fGlobalCT.m[GEMS.rxnXPSID.sort.ind$ix],las=2,ylab="N° rxns per model")

#Find core rxns
nModelsXRxn      <- apply(rxnInclMat==1,1,sum)
nModelsXRxn.perc <- nModelsXRxn/ncol(rxnInclMat)
percsSampleWRxns <- c(seq(0.5,0.90,by = 0.05),
                      seq(0.91,0.99,by = 0.01),
                      seq(0.991,1,by = 0.001))
nRxnsCoreVsPan   <- matrix(0,nrow=2,ncol=length(percsSampleWRxns),
                           dimnames=list(c("N°core","N°pan"),
                                         paste0("Core@",percsSampleWRxns)))
for (percSampleWRxns in percsSampleWRxns){
  rxns.core <- HMR3676.rxnTable[HMR3676.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc>percSampleWRxns]),]
  rxns.pan  <- HMR3676.rxnTable[HMR3676.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc<=percSampleWRxns]),]
  nRxnsCoreVsPan[,paste0("Core@",percSampleWRxns)] <- rbind(nrow(rxns.core),
                                            nrow(rxns.pan))}
barplot(nRxnsCoreVsPan,beside = F,las=2,ylab="N° rxns core vs pan")
###Let's settle for 0.95
percSampleWRxns <- 0.95
rxns.core       <- data.frame(HMR3676.rxnTable[HMR3676.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc>percSampleWRxns]),],"Type"="Core")
rxns.pan.fac    <- data.frame(HMR3676.rxnTable[HMR3676.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc<=percSampleWRxns&nModelsXRxn.perc>0]),],"Type"="Pan","Class"="Facultative")
rxns.pan.abs    <- data.frame(HMR3676.rxnTable[HMR3676.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc==0]),],"Type"="Pan","Class"="Absent")

#Fill the gpr matrix of core rxns with n° model that have a gpr
HMR3675.gpr.num      <- as.matrix(HMR3675.gpr)
HMR3675.gpr.num.core <- HMR3675.gpr.num[row.names(HMR3675.gpr.num)%in%
                                          rxns.core[,"RxnID"],]
HMR3675.gpr.num.core.noemptygene <- HMR3675.gpr.num.core[,
                       apply(HMR3675.gpr.num.core,2,sum)>=1]
nModelsXGPR       <- HMR3675.gpr.num.core.noemptygene*0
nModelsXGene      <- apply(geneInclMat==1,1,sum)
for (gene in colnames(HMR3675.gpr.num.core.noemptygene)){
  nModelsXGPR[HMR3675.gpr.num.core.noemptygene[,gene]!=0,gene] <- nModelsXGene[gene]}

#Compute n° models with a gpr expressed for a core rxn
nModelsXGPR.geneMax      <- apply(nModelsXGPR,2,max,na.rm=T)
nModelsXGPR.geneMax.sort <- sort(nModelsXGPR.geneMax)
a<-annotEntrez[match(names(nModelsXGPR.geneMax.sort),annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
barplot(nModelsXGPR.geneMax.sort,names.arg = a,las=2,
        ylab="N° models with a gpr expressed",main="Core rxns",
        cex.names=0.1)
abline(h=dim(rxnInclMat)*percSampleWRxns,col="red")

#Divide core rxns in univ or context
nModelsXGPR.rxnMaxInAnyGPR <- apply(nModelsXGPR,1,max)
rxns.core.conn <- data.frame(rxns.core[nModelsXGPR.rxnMaxInAnyGPR==0,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR==0],
                             "Class"="Connectivity")
rxns.core.univ <- data.frame(rxns.core[nModelsXGPR.rxnMaxInAnyGPR >= dim(rxnInclMat)[2]*percSampleWRxns,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR >= dim(rxnInclMat)[2]*percSampleWRxns],
                             "Class"="Universal")
rxns.core.cont <- data.frame(rxns.core[nModelsXGPR.rxnMaxInAnyGPR > 0 & nModelsXGPR.rxnMaxInAnyGPR <  dim(rxnInclMat)[2]*percSampleWRxns,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR > 0 & nModelsXGPR.rxnMaxInAnyGPR <  dim(rxnInclMat)[2]*percSampleWRxns],
                             "Class"="Contextual")
write.table(rxns.core.conn,"rxns_core_conn.txt",sep = "\t",quote = F,row.names = F)
write.table(rxns.core.univ,"rxns_core_univ.txt",sep = "\t",quote = F,row.names =F)
write.table(rxns.core.cont,"rxns_core_cont.txt",sep = "\t",quote = F,row.names = F)
write.table(rxns.pan.fac,"rxns_pan_fac.txt",sep = "\t",quote = F,row.names = F)
write.table(rxns.pan.abs,"rxns_pan_abs.txt",sep = "\t",quote = F,row.names = F)

#Check for gpr expression in core rxns 
### Normalize gene inclusion matrix with read count matrix
E.ensg            <- E
colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat))
geneInclMat.c     <- geneInclMat[genes.ensg.common,]
E.ensg.c          <- E.ensg[genes.ensg.common,]
E.ensg.c.zerneg   <- E.ensg.c
E.ensg.c.zerneg[E.ensg.c<0] <- 0
cols <- rainbow(5)

###Example: core rxn - context
nPats.ex           <- 917
rxn.ex             <- "HMR_8604"
rxn.ex.gpr.nModels <- nModelsXGPR[rxn.ex,nModelsXGPR[rxn.ex,]!=0]
rxn.ex.genes       <- colnames(nModelsXGPR)[which(nModelsXGPR[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c[rxn.ex.genes[1],1:nPats.ex],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR3676.rxnTable[HMR3676.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

###Example: core rxn - universal
nPats.ex           <- 917
rxn.ex             <- "HMR_4052"
rxn.ex.gpr.nModels <- nModelsXGPR[rxn.ex,nModelsXGPR[rxn.ex,]!=0]
rxn.ex.genes       <- colnames(nModelsXGPR)[which(nModelsXGPR[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c[rxn.ex.genes[1],1:nPats.ex],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR3676.rxnTable[HMR3676.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

###Example: pan rxn - facultative
nPats.ex           <- 917
rxn.ex             <- "HMR_4429"
rxn.ex.genes       <- colnames(HMR3675.gpr.num)[which(HMR3675.gpr.num[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c[rxn.ex.genes[1],1:nPats.ex],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR3676.rxnTable[HMR3676.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

#Describe core rxns
rxns.subs.min <- 50
rxns.subs <- table(HMR3676.rxnTable[,"Subs"])
rxns.core.subs <- table(rxns.core[,"Subs"])/rxns.subs
rxns.core.conn.subs <- table(rxns.core.conn[,"Subs"])/rxns.subs
rxns.core.univ.subs <- table(rxns.core.univ[,"Subs"])/rxns.subs
rxns.core.cont.subs <- table(rxns.core.cont[,"Subs"])/rxns.subs
rxns.pan.fac.subs   <- table(rxns.pan.fac[,"Subs"])/rxns.subs
rxns.pan.abs.subs   <- table(rxns.pan.abs[,"Subs"])/rxns.subs
rxns.all.subs <-data.frame(rxns.core.conn.subs,
                           rxns.core.cont.subs,
                           rxns.core.univ.subs,
                           rxns.pan.fac.subs,
                           rxns.pan.abs.subs)
rxns.all.subs.f              <- rxns.all.subs[rxns.subs>rxns.subs.min,]
rxns.all.subs.sort           <- rxns.all.subs.f[order(rxns.core.subs[rxns.subs>rxns.subs.min],decreasing = T),]
rownames(rxns.all.subs.sort) <- rxns.all.subs.sort[,1]
rxns.all.subs.sort <- rxns.all.subs.sort[,-grep("Var1",colnames(rxns.all.subs.sort))]
cols <- c("Red","Purple","Dark Red","Blue","Grey")
barplot(as.matrix(t(rxns.all.subs.sort)),cex.names = 0.5,
        beside = F,las=2,horiz=T,col=cols)