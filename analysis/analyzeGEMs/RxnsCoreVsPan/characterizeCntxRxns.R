###Characterize contextual rxns
rm(list=ls())
setwd("~/Documents/Academia/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan")
library(brglm)
require(ConsensusClusterPlus)
require(lmtest)
library(gplots)
library(ggplot2)
library(RColorBrewer)
require(reshape2)
require(varSelRF)
#Param
sign.lev        <- 0.01
min.FC          <- log2(0.5)
percSampleWRxns <- 0.95
sign.ct         <- 0.001
minFSampleXMut  <- 0.10

#Load data
load("../../data/annotation.RData")
load("../../data/preprocessed.RData") 
  rxnIncMatFile     <- c("../rxnInclMat.txt")
  rxnInclMat        <- read.delim(rxnIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat  <- rxnInclMat[,-grep("X.1",colnames(rxnInclMat))]
  geneIncMatFile    <- c("../geneInclMat.txt")
  geneInclMat       <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
geneInclMat <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
  HMR2.gpr          <- read.delim("../../docs/cHMR3765_rxnGeneMatrix.txt")
  rownames(HMR2.gpr)<-HMR2.gpr[,1]
HMR2.gpr      <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]
HMR2.rxnTable <- read.delim("../../data/cHMR3765_rxnTable.txt",stringsAsFactors=F)
cntx.rxns     <- read.delim("rxns_cntx.txt",stringsAsFactors=F)

#1: Preprocess
PSID.m        <- colnames(rxnInclMat)
  lib.size          <- with(y$samples, lib.size * norm.factors)
  E                 <- t(log2(t(y$counts + 0.5)/(lib.size + 1) * 1e+06))
  E.ensg            <- E
  colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
  rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
  genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat))
geneInclMat.c <- geneInclMat[genes.ensg.common,]
E.ensg.c      <- E.ensg[genes.ensg.common,] #n common genes x nS (1082)
HMR2.gpr.num  <- as.matrix(HMR2.gpr)

#2: Create a factor
  PSID         <- rownames(y$samples)
  PSID.s       <- gsub(pattern="-",replacement=".",x=PSID)
  filter       <- PSID.s%in%colnames(rxnInclMat)
  mapper       <- match(colnames(rxnInclMat),PSID.s[filter])
fGlobalCT.m  <- fGlobalCT[filter][mapper]
  names(fGlobalCT.m)  <- PSID.s[filter][mapper]
  fMuts.f      <- fMuts[filter,]
  fMuts.m      <- fMuts.f[mapper,]
mut          <- data.matrix(data.frame(fMuts.m))-1
mut.freq     <- mut[,apply(mut,2,sum)>minFSampleXMut*nrow(mut)]
E.ensg.c.m   <- E.ensg.c[,filter][,mapper] #n common genes x nS in models (917)
cols         <- colorRampPalette(brewer.pal(12, "Set3"))(nlevels(fGlobalCT.m))
  names(cols)  <- levels(fGlobalCT.m)
colScale     <- scale_fill_manual(values = cols)

#3: Example for gpr expression in core rxns 
### Normalize gene inclusion matrix with read count matrix
E.ensg.c.zerneg             <- E.ensg.c
E.ensg.c.zerneg[E.ensg.c<0] <- 0
cols5 <- brewer.pal(5,"Set2")

###Example: cntx rxn
rxn.ex             <- "HMR_1037"
rxn.ex.genes       <- colnames(HMR2.gpr.num)[which(HMR2.gpr.num[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c.zerneg[rxn.ex.genes[1],],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols5,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR2.rxnTable[HMR2.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

#4: How are the cntx rxns scattered across models?
# rxnInclMat.cntx   <- rxnInclMat[match(rownames(rxnInclMat),cntx.rxns$Eq)
#                                 [!is.na(match(rownames(rxnInclMat),cntx.rxns$Eq))],]
# hm  <- heatmap.2(as.matrix(rxnInclMat.cntx),Rowv = T,Colv = T,symm = F,
#                  scale="none",col=c("White","Black"),density.info="none",
#                  trace="none",labCol=NA,labRow=NA)

#4: Test if presence of a cntx in a model is attributable to a certain factor ~10'
if (!file.exists("rxns_cntx_inclusiontestres.RData")){
  p.incl     <- matrix(nrow=nrow(cntx.rxns),ncol=5,
                   dimnames=list(cntx.rxns$RxnID,c("pMut","pCT","AICmut","AICct","AICnull")))
  k <- 0
  for (rxn in cntx.rxns$RxnID){
    k <- k+1
    print(paste0("Doing ",k,"/",length(cntx.rxns$RxnID)))
    rxn.eq           <- HMR2.rxnTable[HMR2.rxnTable$RxnID==rxn,"Eq"]
    ### Test with a linear model if inclusion goes with ct or mut
    fit.mut          <- brglm(as.numeric(rxnInclMat[rxn.eq,])~fGlobalCT.m+mut.freq)
    fit.ct           <- update(fit.mut, ~.- mut.freq)
    test.mut         <- lrtest(fit.ct,fit.mut)
    fit.null         <- update(fit.ct, ~.- fGlobalCT.m)
    test.ct          <- lrtest(fit.null,fit)
    p.incl[rxn,]     <- c(test.mut$"Pr(>Chisq)"[2],test.ct$"Pr(>Chisq)"[2],
                          fit.mut$aic,fit.ct$aic,fit.null$aic)}} else {
  load("rxns_cntx_inclusiontestres.RData")}

###Sort 
p      <- p.incl[,c("pMut","pCT")]
p.adj  <- apply(p,2,p.adjust,"fdr")
p.sort <- p.adj[order(p.adj[,"pCT"]),]

###Compute for interesting rxns the fraction in which a ct is represented (where p < 0.01)
cntx.rxns.int  <- cntx.rxns[p.adj[,"pCT"] < sign.ct,]
cntx.rxns.nxct <- matrix(0,nrow=nlevels(fGlobalCT.m),ncol=dim(cntx.rxns.int)[1],
                            dimnames=list(levels(fGlobalCT.m),cntx.rxns.int$Eq))
k <- 0
for (rxn in cntx.rxns.int$Eq){
  k <- k+1
  print(paste0("Doing: #",k))
  a <- t(rxnInclMat[rxn,])
  rxn.inter.sumXFactor <- aggregate(a~fGlobalCT.m,FUN = sum,simplify = T)
  cntx.rxns.nxct[,rxn] <- rxn.inter.sumXFactor[,rxn]
#   barplot(rxn.inter.sumXFactor[,2]/table(fGlobalCT.m),names.arg = rxn.inter.sumXFactor[,1] ,
#           las=2,main=rxn,ylim = c(0,1))
}
cntx.rxns.fracxct <- t(sweep(cntx.rxns.nxct,1,table(fGlobalCT.m),"/"))

###Consensus clustering on the partitioning of a reaction among ct
ccp     <- ConsensusClusterPlus(cntx.rxns.fracxct,maxK=5,reps=1000,pItem=0.8,pFeature=1,
                            title="Fraction of models with the reaction",clusterAlg="hc",
                            distance="pearson",
                            seed=1262118388.71279,plot=NULL)
clust4k <- ccp[[4]]
cntx.rxns.fracxct.clust4k <- cntx.rxns.fracxct[,clust4k$consensusTree$order]
hmcols  <- colorRampPalette(brewer.pal(9, "GnBu"))(256)
hm      <- heatmap.2(cntx.rxns.fracxct.clust4k,Rowv = T,Colv = NA,symm = F,
                    scale="none",col=hmcols,density.info="none",trace="none",
                    labRow=NA,dendrogram="row",labCol=toupper(colnames(cntx.rxns.fracxct.clust4k)),
                    ColSideColors=brewer.pal(4, "Set3")
                    [clust4k$consensusClass[clust4k$consensusTree$order]])

#Order reaction by maximal variability
var.score <- rank(apply(cntx.rxns.fracxct,1,mad))
p.score   <- rank(p.adj[p.adj[,"pCT"]<sign.ct,"pCT"])
tot.score <- 0.2*p.score+0.8*var.score
cntx.rxns.fracxct.top <- data.frame(cntx.rxns.int
                                    [order(tot.score),],
                                    score=tot.score[order(tot.score)])
for (rxn in cntx.rxns.fracxct.top$Eq[1:2]){
  a <- t(rxnInclMat[rxn,])
  rxn.inter.sumXFactor <- aggregate(a~fGlobalCT.m,FUN = sum,simplify = T)
  barplot(rxn.inter.sumXFactor[,2]/table(fGlobalCT.m),names.arg = rxn.inter.sumXFactor[,1] ,
      las=2,main=rxn,ylim = c(0,1),col=cols)
}

###Select top 5 rxns per cluster
if (!file.exists("rxns_cntx_varsel.RData")){
  features <- t(rxnInclMat[match(cntx.rxns.fracxct.top[1:2000,]$Eq,rownames(rxnInclMat))
                           [!is.na(match(cntx.rxns.fracxct.top[1:2000,]$Eq,rownames(rxnInclMat)))],])
    colnames(features) <- cntx.rxns.fracxct.top[1:2000,]$RxnID
  classes  <- as.character(fGlobalCT.m)
  for (ct in levels(fGlobalCT.m)){
    classes[which(classes==ct)] <- clust4k$consensusClass[ct]}
    classes <- as.factor(classes)
  var.sel  <- varSelRF(features,classes)
  save("var.sel",file="rxns_cntx_varsel.RData")} else {load("rxns_cntx_varsel.RData")}
top.rxns.eq                      <- HMR2.rxnTable[HMR2.rxnTable$RxnID %in% var.sel$selected.vars,"Eq"]
cntx.rxns.fracxct.clust4k.varsel <- cntx.rxns.fracxct.clust4k[top.rxns.eq,]

###Filter out identical reactions
filter  <- duplicated(apply(cntx.rxns.fracxct.clust4k.varsel,1,sum))
cntx.rxns.fracxct.clust4k.varsel.f <- cntx.rxns.fracxct.clust4k.varsel[!filter,]
hm      <- heatmap.2(cntx.rxns.fracxct.clust4k.varsel.f,Rowv = T,Colv = NA,symm = F,
                     scale="none",col=hmcols,density.info="none",trace="none",
                     dendrogram="row",
                     ColSideColors=brewer.pal(4, "Set3")
                     [clust4k$consensusClass[clust4k$consensusTree$order]])
tmp        <- rownames(cntx.rxns.fracxct.clust4k.varsel.f[hm$rowInd,])
cntx_clust <- HMR2.rxnTable[match(tmp,HMR2.rxnTable$Eq),"RxnID"]

###Plot top rxns most important
top.rxns.examples <- c("HMR_0450","HMR_9611","HMR_4187","HMR_4169","HMR_4160",
                       "HMR_4843","HMR_6713","HMR_4101","HMR_4296","HMR_1037")
for (rxn in top.rxns.examples){
  ### what genes are associated?
  rxn.genes   <- colnames(HMR2.gpr)[which(HMR2.gpr[rxn,]!=0)]
  rxn.genes   <- rxn.genes[rxn.genes %in% genes.ensg.common]
  rxn.genes.n <- length(rxn.genes)
  E.rxn.genes           <- as.matrix(E.ensg.c.m[rxn.genes,])
  if (dim(E.rxn.genes)[1]<dim(E.rxn.genes)[2]){E.rxn.genes <- t(E.rxn.genes)}
  colnames(E.rxn.genes) <- annotEntrez[match(rxn.genes,
                                           annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
  fGlobalCT.m.o         <- factor(fGlobalCT.m,
                                  levels = levels(fGlobalCT.m)
                                            [clust4k$consensusTree$order])
  E.rxn.genes.annot   <- data.frame(E.rxn.genes,fGlobalCT.m.o)
  dfmelt <- melt(E.rxn.genes.annot, measure.vars = 1:rxn.genes.n)

  filename <- paste0("Plots/",rxn,"_boxplotCT.pdf")
  p1 <- ggplot(dfmelt, aes(x=fGlobalCT.m.o, y=value, fill=fGlobalCT.m.o))+
    geom_boxplot()+
    facet_grid(.~variable)+
    xlab("")+
    ylab("")+
    guides(fill=F)+
    theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))+
    ggtitle(paste(cntx.rxns[match(rxn,cntx.rxns$RxnID),"Eq"],
                  collapse="\n"))+
    theme(plot.title = element_text(size = rel(1)))
  p1 <- p1 + colScale
  ggsave(filename,plot = p1,width = 13.8,height = 8.62)}

###Compute for interesting rxns the fraction in which a mut is represented (where p < 0.01)
cntx.rxns.int.mut  <- cntx.rxns[p.adj[,"pMut"] < sign.ct,]
rxn.eq             <- cntx.rxns.int.mut$Eq
fit.mut            <- brglm(as.numeric(rxnInclMat[rxn.eq,])~fGlobalCT.m+mut.freq)
rxn.mut.factor     <- factor(mut.freq[,"TP53"])
rxn.inter.sumXFactor <- aggregate(t(rxnInclMat[rxn.eq,])~rxn.mut.factor,FUN = sum,simplify = T)
pdf("Plots/HMR_9802_boxplotTP53.pdf")
barplot(rxn.inter.sumXFactor[,2]/table(rxn.mut.factor),
        names.arg = paste("TP53",rxn.inter.sumXFactor[,1]),
        main=rxn.eq,ylim = c(0,1),col=c("white","dark grey"))
  rxn.genes   <- colnames(HMR2.gpr)[which(HMR2.gpr[cntx.rxns.int.mut$RxnID,]!=0)]
  rxn.genes   <- rxn.genes[rxn.genes %in% genes.ensg.common]
  E.rxn.genes <- as.matrix(E.ensg.c.m[rxn.genes,])
    if (dim(E.rxn.genes)[1]<dim(E.rxn.genes)[2]){E.rxn.genes <- t(E.rxn.genes)}
    colnames(E.rxn.genes) <- annotEntrez[match(rxn.genes,
                                             annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
  E.rxn.genes.annot <- data.frame(E.rxn.genes,mut=rxn.mut.factor)
  dfmelt            <- melt(E.rxn.genes.annot, measure.vars = 1:rxn.genes.n)
  p1 <- ggplot(dfmelt, aes(x=mut, y=value,fill=mut))+
    geom_boxplot()+
    facet_grid(.~variable)+
    xlab("TP53")+
    ylab("")+
    guides(fill=F)+
    theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))+
    ggtitle(rxn.eq)+
    theme(plot.title = element_text(size = rel(1)))
  p1
dev.off()

#7: Save
writeLines(cntx_clust,"rxns_cntx_clust.txt")

#8: Test if presence of a rxn correlates with mutation burden
mut.burden <- apply(mut,1,sum)
p.incl     <- matrix(nrow=nrow(cntx.rxns),ncol=2,
                     dimnames=list(cntx.rxns$RxnID,c("pMutBurd",
                                                     "coeffMutBurd")))
k <- 0
for (rxn in cntx.rxns$RxnID){
  k <- k+1
  print(paste0("Doing ",k,"/",length(cntx.rxns$RxnID)))
  rxn.eq           <- HMR2.rxnTable[HMR2.rxnTable$RxnID==rxn,"Eq"]
  ### Test with a linear model if inclusion goes with ct or mut
  fit.mut          <- brglm(as.numeric(rxnInclMat[rxn.eq,])~0+
                              fGlobalCT.m+mut.burden)
  fit.ct           <- update(fit.mut, ~.- mut.burden)
  test.mut         <- lrtest(fit.ct,fit.mut)
  p.incl[rxn,]     <- c(test.mut$"Pr(>Chisq)"[2],
                        fit.mut$coefficients["mut.burden"])}

###Sort 
p.adj    <- p.adjust(p.incl[,"pMutBurd"],"fdr")
p.sort   <- sort(p.adj)

#Select rxns whose loss or gain is significantly associated with mutation burden
rxn.sign.loss   <- cntx.rxns$RxnID[p.adj < sign.lev & p.incl[,"coeffMutBurd"]<0]
rxn.sign.gain   <- cntx.rxns$RxnID[p.adj < sign.lev & p.incl[,"coeffMutBurd"]>0]

#There are only losses. Plot density per mutation burdern of the reaction: equally absent/present?
rxn.sign.loss.eq <- HMR2.rxnTable[HMR2.rxnTable$RxnID%in%rxn.sign.loss,"Eq"]
rxnInclMat.loss <- rxnInclMat[rownames(rxnInclMat) %in% rxn.sign.loss.eq,]
rxnInclMat.loss.mb <- c()
for(rxn in rownames(rxnInclMat.loss)){
  rxnInclMat.loss.mb   <- rbind(rxnInclMat.loss.mb,
                                data.frame(factor(ifelse(rxnInclMat.loss[rxn,]==1,
                                                 "Present","Absent")),
                                   mut.burden))}
  colnames(rxnInclMat.loss.mb) <- c("Rxn.status","Mut.burden")
rxnInclMat.loss.melt <- melt(rxnInclMat.loss.mb,
                             measure.vars = "Mut.burden")
ggplot(rxnInclMat.loss.melt, aes(x=value, fill=Rxn.status)) + 
  geom_density(alpha=.3) +
  xlab("Mutation burden")

#Provide an example
rxnInclMat.loss.mb.ex1 <- data.frame(factor(ifelse(rxnInclMat.loss[1,]==1,
                         "Present","Absent")),mut.burden)
colnames(rxnInclMat.loss.mb.ex1) <- c("Rxn.status","Mut.burden")
rxnInclMat.loss.melt.ex1 <- melt(rxnInclMat.loss.mb.ex1,
                             measure.vars = "Mut.burden")
ggplot(rxnInclMat.loss.melt.ex1, aes(x=value, fill=Rxn.status)) + 
  geom_density(alpha=.3) +
  xlab("Mutation burden") +
  ggtitle(rxn.sign.loss.eq[1])

#Any evidence at the expression level? No
rxn.sign.loss.genes <- c()
for (rxn in rxn.sign.loss){
  ### what genes are associated?
  rxn.genes   <- colnames(HMR2.gpr)[which(HMR2.gpr[rxn,]!=0)]
  rxn.genes   <- rxn.genes[rxn.genes %in% genes.ensg.common]
  rxn.sign.loss.genes <- c(rxn.sign.loss.genes,rxn.genes)}
rxn.sign.loss.genes <- unique(rxn.sign.loss.genes)
writeLines(rxn.sign.loss.genes,"genesasswrxnlosswmutburden_ENSG.txt")

#Is this spurious? Is it robust to permutation test of the mut burden?

#Is the mut burden corr with the n° of rxns? Yes, but very limited effect
nrxnsxmodel <- apply(rxnInclMat==1,2,sum)
plot(log2(mut.burden),nrxnsxmodel)
cortest <- cor.test(mut.burden,nrxnsxmodel)
fit     <- glm(nrxnsxmodel~fGlobalCT.m+mut.burden+fGlobalCT.m:mut.burden)
