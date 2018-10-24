###Compute similarity between models and cluster
rm(list=ls())
setwd("~/Box Sync/3rd Semester (May 2013 - November 2013)/Projects/RxnLands/Codes/analyzeGEMs/ModelSimilarity")
require(ade4)
require(ConsensusClusterPlus)
require(doMC)
require(glmnet)
library(gplots)
require(randomForest)
require(RColorBrewer)
library(ROCR)
library(stringdist)
source("pcaElInclMatrix.R")

#Load and params
load("../../data/annotation.RData")
load("../../data/preprocessed.RData")
HMR2.rxnTable  <- read.delim("../../docs/cHMR3765_rxnTable.txt")
  HMR2.gpr           <- read.delim("../../docs/cHMR3765_rxnGeneMatrix.txt")
  rownames(HMR2.gpr) <-HMR2.gpr[,1]
HMR2.gpr       <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]
elIncMatFiles  <- c("../rxnInclMat.txt")#,"../geneInclMat.txt","../metInclMat.txt")
geneIncMatFile <- c("../geneInclMat.txt")#,"../geneInclMat.txt","../metInclMat.txt")
  geneInclMat        <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
  geneInclMat        <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
registerDoMC(cores=4)

#Compute hamming distance for each element incl mat: focus on rxns
H         <- list()
elInclMat <- list()
for (elIncMatFile in elIncMatFiles){
  el <- strsplit(strsplit(elIncMatFile,"/")[[1]][2],"I")[[1]][1]
  elInclMat <- read.delim(elIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
  elInclMat[[el]] <- elInclMat[,-grep("X.1",colnames(elInclMat))]
  nEls      <- dim(elInclMat[[el]])[1]
  H.notNorm <- stringdistmatrix(elInclMat[[el]],elInclMat[[el]],method="hamming")
  H[[el]]   <- H.notNorm*(-1/nEls)+1 #Scale between 0 (very dissimilar) to 1 (identical)
  dimnames(H[[el]])<-list(colnames(elInclMat[[el]]),colnames(elInclMat[[el]]))
}
rxnInclMat   <- elInclMat[["rxn"]]

#Create a factor
PSID         <- rownames(y$samples)
PSID.m       <- colnames(rxnInclMat)
PSID.s       <- gsub(pattern="-",replacement=".",x=PSID)
filter       <- PSID.s%in%PSID.m
mapper       <- match(PSID.m,PSID.s[filter])
fGlobalCT.m  <- fGlobalCT[filter][mapper]
names(fGlobalCT.m)  <- PSID.s[filter][mapper]
nCancerTypes <-nlevels(fGlobalCT.m)

#PCA
par(mfrow=c(2,1))
pcaElInclMatrix(geneInclMat)
pcaElInclMatrix(rxnInclMat)
par(mfrow=c(1,1))

#Create expression matrix
  lib.size          <- with(y$samples, lib.size * norm.factors)
  E                 <- t(log2(t(y$counts + 0.5)/(lib.size + 1) * 1e+06))
  E.ensg            <- E
  colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
  rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
  genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat))
geneInclMat.c     <- geneInclMat[genes.ensg.common,]
E.ensg.c          <- E.ensg[genes.ensg.common,PSID.m]

#Consensus clustering upon distance
if (!file.exists("consensusClustRxnSimil_reps1000.RData")){
ccp <- ConsensusClusterPlus(H[["rxn"]],maxK=5,reps=1000,pItem=0.8,pFeature=1,
                               title="rxnSimilarity_consensus_cluster",clusterAlg="hc",
                               distance="pearson",
                               seed=1262118388.71279,plot=NULL)
icl <- calcICL(results,title="rxnSimilarity_consensus_cluster",
               plot=NULL,writeTable=FALSE)} else {load("consensusClustRxnSimil_reps1000.RData")}
clust3k     <- ccp[[3]]
H.clust3k   <- H[["rxn"]][clust3k$consensusTree$order,
                        clust3k$consensusTree$order]
hmcols <- colorRampPalette(brewer.pal(3, "YlGnBu"))(256)
hm     <- heatmap.2(H.clust3k,Rowv = NA,Colv = NA,symm = F,scale="none",col=hmcols,density.info="none",trace="none",labCol=NA,labRow=NA)

#Is belonging to class 1 or 3 (response) due to some rxn inclusions (obs)?
model.class          <- clust3k$consensusClass
rxnInclMat.class1or3 <- t(rxnInclMat[,model.class!=2])
class1or3            <- as.factor(model.class[model.class!=2])
names(class1or3)     <- names(model.class[model.class!=2])

###Lasso selection of coefficient most predictive of class belonging
out.cv      <- cv.glmnet(rxnInclMat.class1or3,class1or3,parallel = T,
                 family = "binomial", nfolds=20)
coef.min    <- coef(out.cv, s = "lambda.1se")
active.min  <- which(coef.min != 0)
index.min   <- coef.min[active.min]
rxns.eq.sel <- rownames(coef.min)[active.min][-1] #Remove "Intercept"
rxnInclMat.class1or3.sel <- t(rxnInclMat[rxns.eq.sel,model.class!=2])

###Test if there is any significant difference in the two classes in the inclusion of the lasso rxns
wilcox      <- matrix(0,ncol=1,nrow=ncol(rxnInclMat.class1or3.sel),
                      dimnames=list(colnames(rxnInclMat.class1or3.sel),"WilcoxP"))
for (rxn in colnames(rxnInclMat.class1or3.sel)){
  res <- wilcox.test(as.numeric(rxnInclMat.class1or3.sel[,rxn])~class1or3)
  wilcox[rxn,] <- res$p.value
}
wilcox[,1] <- p.adjust(wilcox,method = "fdr")

###Retrieve significantly different rxns and cluster
rxns.eq.sel2 <- rownames(wilcox)[wilcox<0.01]
rxnInclMat.clust3k.sel <- rxnInclMat[rxns.eq.sel2,
                                     clust3k$consensusTree$order]
hm          <- heatmap.2(as.matrix(rxnInclMat.clust3k.sel),
                         Rowv = NA,Colv = F,symm = F,scale="none",
                         col=c("White","Black"),density.info="none",
                         trace="none",labCol=NA)

#What are the genes beyond the rxns that make class 1 and 3 different?
rxns.table.sel2  <- HMR2.rxnTable[match(rxns.eq.sel2,HMR2.rxnTable$Eq),]
genes.sel2       <- colnames(HMR2.gpr)[which(apply(
                    HMR2.gpr[match(rxns.table.sel2$RxnID,rownames(HMR2.gpr))
                             ,],2,sum)>0)]
PSID.class1or3   <- names(class1or3)
geneInclMat.sel2 <- geneInclMat.c[genes.sel2[genes.sel2%in%rownames(geneInclMat.c)],]
hm <- heatmap.2(as.matrix(geneInclMat.sel2[,clust3k$consensusTree$order]),
                trace="none",
                scale = "none",col=c("White","Black"),Rowv=T,Colv=F,symm=F,
                labCol=NA,labRow=annotEntrez[genes.sel2.entz,"hgnc_symbol"],
                ColSideColors=as.character(model.class[clust3k$consensusTree$order]))

E.sel2 <- E.ensg.c[genes.sel2[genes.sel2%in%rownames(E.ensg.c)],]
hm <- heatmap.2(E.sel2[,clust3k$consensusTree$order],trace="none",col=hmcols,
                Rowv=T,Colv=F,symm=F,scale = "none",symbreaks=F,
                breaks=seq(-2,10,12/length(hmcols)),
                labCol=NA,labRow=annotEntrez[genes.sel2.entz,"hgnc_symbol"],
                ColSideColors=as.character(model.class[clust3k$consensusTree$order]))

#Are there genes NOT encoding for these rxns relevant?!
y.class1or3 <- y[,rownames(y$samples)%in%PSID.class1or3]
v.fit       <- voom(y.class1or3,model.matrix(~class1or3),plot=T)
lm.fit      <- lmFit(v.fit,model.matrix(~class1or3))
eB          <- eBayes(lm.fit)
tT          <- topTable(eB,number = Inf)
genes.DE.class1or3.ensg <- tT$ensembl_gene_id[tT$adj.P.Val < 0.05]
genes.DE.class1or3.hgnc <- tT$hgnc_symbol[tT$adj.P.Val < 0.05]
E.sel3      <- E.ensg.c[genes.DE.class1or3.ensg
                        [genes.DE.class1or3.ensg%in%rownames(E.ensg.c)],]
hm <- heatmap.2(E.sel3[,clust3k$consensusTree$order],trace="none",col=hmcols,
                Rowv=T,Colv=F,symm=F,scale = "none",symbreaks=F,
                labCol=NA,labRow=genes.DE.class1or3.hgnc[genes.DE.class1or3.ensg%in%rownames(E.ensg.c)],
                ColSideColors=as.character(model.class[clust3k$consensusTree$order]))

