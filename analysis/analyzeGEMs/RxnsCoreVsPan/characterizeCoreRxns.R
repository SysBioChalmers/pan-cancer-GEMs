###Characterize core contextuale rxns
rm(list=ls())
setwd("~/Box Sync/3rd Semester (May 2013 - November 2013)/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan")
library(brglm)
require(lmtest)
library(gplots)
library(ggplot2)
library(RColorBrewer)

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
HMR2.gpr    <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]
core.rxns   <- read.delim("rxns_core.txt",stringsAsFactors=F)

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

#3: Distinguish pan iso conn by counting the extent to which a gpr is shared
#   (95%: pan ~ <95%: iso ~ 0%: conn)
### Fill the gpr matrix of core rxns with n° model that have a gpr
HMR2.gpr.num      <- as.matrix(HMR2.gpr)
HMR2.gpr.num.core <- HMR2.gpr.num[row.names(HMR2.gpr.num)%in%
                                    core.rxns[,"RxnID"],]
HMR2.gpr.num.core.noemptygene <- HMR2.gpr.num.core[, #This excludes the conn
                                        apply(HMR2.gpr.num.core,2,sum)>=1]
nModelsXGPR       <- HMR2.gpr.num.core.noemptygene*0
nModelsXGene      <- apply(geneInclMat==1,1,sum)
for (gene in colnames(HMR2.gpr.num.core.noemptygene)){
  nModelsXGPR[HMR2.gpr.num.core.noemptygene[,gene]!=0,gene] <- nModelsXGene[gene]}

### Compute n° models with a gpr expressed for a core rxn
nModelsXGPR.geneMax      <- apply(nModelsXGPR,2,max,na.rm=T)
nModelsXGPR.geneMax.sort <- sort(nModelsXGPR.geneMax)
label <-annotEntrez[match(names(nModelsXGPR.geneMax.sort),
                     annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
barplot(nModelsXGPR.geneMax.sort,names.arg = label,las=2,
        ylab="N° models with a gpr expressed",main="Core rxns",
        cex.names=0.1)
abline(h=dim(rxnInclMat)*percSampleWRxns,col="red")

### Divide core rxns in pan, iso, conn
nModelsXGPR.rxnMaxInAnyGPR <- apply(nModelsXGPR,1,max)
threshold      <- dim(rxnInclMat)[2]*percSampleWRxns
rxns.core.conn <- data.frame(core.rxns[nModelsXGPR.rxnMaxInAnyGPR==0,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR==0],
                             "Class"="Connectivity")
rxns.core.pan  <- data.frame(core.rxns[nModelsXGPR.rxnMaxInAnyGPR >= threshold,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR >= dim(rxnInclMat)[2]*percSampleWRxns],
                             "Class"="Pan")
rxns.core.iso  <- data.frame(core.rxns[nModelsXGPR.rxnMaxInAnyGPR > 0 & nModelsXGPR.rxnMaxInAnyGPR <  threshold,],
                             "nModelsXGPR.rxnMaxInAnyGPR"=nModelsXGPR.rxnMaxInAnyGPR[nModelsXGPR.rxnMaxInAnyGPR > 0 & nModelsXGPR.rxnMaxInAnyGPR <  dim(rxnInclMat)[2]*percSampleWRxns],
                             "Class"="Iso")
#4: Bootstrap statistics
nbootstraps <- 1000
n.pan     <- function(modelsxgpr,i) sum(modelsxgpr[i]>threshold)
n.pan.bs  <- boot(nModelsXGPR.rxnMaxInAnyGPR,statistic = n.pan,R = nbootstraps)
n.pan.ci  <- boot.ci(n.pan.bs,type="basic")
n.iso     <- function(modelsxgpr,i) sum(modelsxgpr[i]<threshold & modelsxgpr[i]>0)
n.iso.bs  <- boot(nModelsXGPR.rxnMaxInAnyGPR,statistic = n.iso,R = nbootstraps)
n.iso.ci  <- boot.ci(n.iso.bs,type="basic")
n.conn    <- function(modelsxgpr,i) sum(modelsxgpr[i]==0)
n.conn.bs <- boot(nModelsXGPR.rxnMaxInAnyGPR,statistic = n.conn,R = nbootstraps)
n.conn.ci <- boot.ci(n.conn.bs,type="basic")

n   <- c(n.pan.bs$t0,n.iso.bs$t0,n.conn.bs$t0)
n.l <- c(n.pan.ci$basic[1,4],n.iso.ci$basic[1,4],n.conn.ci$basic[1,4])
n.u <- c(n.pan.ci$basic[1,5],n.iso.ci$basic[1,5],n.conn.ci$basic[1,5])

barplot2(n,horiz=F,beside=F,las=2,names.arg = c("Pan","Iso","Connectivity"),         
         ci.l =n.l,ci.u =n.u,
         col = brewer.pal(length(n),"YlOrRd")[c(3,2,1)],
         ylab="N° reactions",plot.ci=T,
         legend.text=F,ylim=c(0,4000))

#5: Give examples for gpr expression in core rxns 
E.ensg.c.zerneg               <- E.ensg.c.m
E.ensg.c.zerneg[E.ensg.c.m<0] <- 0
cols5 <- brewer.pal(5,"Set2")

###Example: core rxn - iso
rxn.ex             <- "HMR_8604"
rxn.ex.gpr.nModels <- nModelsXGPR[rxn.ex,nModelsXGPR[rxn.ex,]!=0]
rxn.ex.genes       <- colnames(nModelsXGPR)[which(nModelsXGPR[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c.m[rxn.ex.genes[1],],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols5,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR2.rxnTable[HMR2.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

###Example: core rxn - universal
rxn.ex             <- "HMR_4052"
rxn.ex.gpr.nModels <- nModelsXGPR[rxn.ex,nModelsXGPR[rxn.ex,]!=0]
rxn.ex.genes       <- colnames(nModelsXGPR)[which(nModelsXGPR[rxn.ex,]!=0)]
rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
rxn.ex.gene.sort   <- sort.int(E.ensg.c.m[rxn.ex.genes[1],],index.return = T)

barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
        col=cols5,legend.text = rxn.ex.genes.hgnc,horiz=T,
        main=HMR2.rxnTable[HMR2.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
        xlab="Size-corrected log-cpm")

#6: Describe iso - most abundant gpr in factor per iso rxn
##First cluster "identical" rxns
HMR2.gpr.iso   <- HMR2.gpr.num[rownames(HMR2.gpr.num)%in%rxns.core.iso$RxnID,]
HMR2.gpr.iso.f <- HMR2.gpr.iso[,apply(HMR2.gpr.iso!=0,2,any)]
hc.rows        <- hclust(dist(HMR2.gpr.iso.f))
rxn.clusters   <- cutree(hc.rows, h=0.5) #Below 1 all rxns in clusters are "identical" (same gene annotation)
n.clusters     <- max(rxn.clusters)
###Then check for each rxn cluster if any gene associated is ct or mut dependent
p                      <- matrix(nrow=1,ncol=4,
                         dimnames=list("",c("rxn","gene","Chi^2p_mut","Chi^2p_ct")))
rxn.maxgprxct          <- matrix(0,nrow=n.clusters,ncol=nlevels(fGlobalCT.m),
                         dimnames=list(NULL,levels(fGlobalCT.m)))
rxn.maxgprxct.cellnote <- matrix("",nrow=n.clusters,ncol=nlevels(fGlobalCT.m),
                         dimnames=list(NULL,levels(fGlobalCT.m)))
rxn.repr.clusters      <- matrix("",nrow=n.clusters,ncol=1)
for (rxn.serial in 1:n.clusters){
  rxns.incluster <- names(which(rxn.clusters==rxn.serial))
  #Pick the first (random) rxn in the cluster
  rxn            <- rxns.incluster[1]
  rxn.repr.clusters[rxn.serial] <- rxn
  print(paste0("Doing",rxn))
  ### what genes are associated?
  rxn.genes <- colnames(HMR2.gpr)[which(HMR2.gpr[rxn,]!=0)]
  rxn.genes <- rxn.genes[rxn.genes %in% genes.ensg.common]
  ### Does it vary with the cancer type?
  rxn.genes.n <- length(rxn.genes)
  if (rxn.genes.n > 5){
    print(paste0(rxn," has ",rxn.genes.n," genes associated. Skipping to next rxn..."))
    next
  }
  
  ### Test with a linear model if inclusion goes with ct
  p.rxn <- matrix(0,nrow=length(rxn.genes),ncol=2,dimnames=list(rxn.genes,c("pMUT","pCT")))
  for (gene in rxn.genes){
    fit.mut  <- brglm(as.numeric(geneInclMat.c[gene,])~fGlobalCT.m+mut.freq)
    fit.ct   <- update(fit.mut, ~.- mut.freq)
    test.mut <- lrtest(fit.ct,fit.mut)
    fit.null <- update(fit.ct, ~.- fGlobalCT.m)
    test.ct  <- lrtest(fit.null,fit)
    p        <- rbind(p,c(rxn,
                          annotEntrez[match(gene,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"],
                          test.mut$"Pr(>Chisq)"[2],test.ct$"Pr(>Chisq)"[2]))
    p.rxn[gene,] <- c(test.mut$"Pr(>Chisq)"[2],test.ct$"Pr(>Chisq)"[2])
  }
  if (any(p.rxn[,"pMUT"] < sign.ct)){
    print(paste0(rxn," encoded by at least a mut dependent gene!"))}
  if (any(p.rxn[,"pCT"] < sign.ct)){
    print(paste0(rxn," encoded by at least a ct dependent gene!"))
    ### If cT specific, bin data according to cancer type
    E.rxn.genes       <- E.ensg.c.m[rxn.genes,]
    rownames(E.rxn.genes) <- annotEntrez[match(rownames(E.rxn.genes),
                                annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
    fGlobalCT.m.o     <- factor(fGlobalCT.m,levels = levels(fGlobalCT.m)[order(tT[,1:12])])
    E.rxn.genes.annot <- data.frame(t(E.rxn.genes),fGlobalCT.m.o)
    dfmelt <- melt(E.rxn.genes.annot, measure.vars = 1:rxn.genes.n)
  
    filename <- paste0("Plots/",rxn,"_boxplotCT.pdf")
    p1 <- ggplot(dfmelt, aes(x=fGlobalCT.m.o, y=value, fill=fGlobalCT.m.o))+
           geom_boxplot()+
           facet_grid(.~variable)+
           xlab("")+
           ylab("")+
           guides(fill=F)+
           theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))+
           ggtitle(paste(rxns.core.iso[match(rxns.incluster,rxns.core.iso$RxnID),"Eq"],
                         collapse="\n"))+
           theme(plot.title = element_text(size = rel(0.8)))
    p1 <- p1 + colScale
    p2 <- ggplot(dfmelt, aes(x=variable, y=value, fill=variable))+
      geom_boxplot()+
      facet_grid(.~fGlobalCT.m.o)+
      theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))+
      ggtitle(rxns.core.iso[rxn,"Eq"])
    ggsave(filename,plot = p1,width = 13.8,height = 8.62)
    ### Get the most abundant in each ct
    for (ct in levels(fGlobalCT.m)){
      nmodels.ct.xgene <- apply(geneInclMat.c[rxn.genes,fGlobalCT.m==ct],1,sum)
      maxgene.ct.ind   <- which.max(nmodels.ct.xgene)
      rxn.maxgprxct[rxn.serial,ct]          <- maxgene.ct.ind
      rxn.maxgprxct.cellnote[rxn.serial,ct] <- annotEntrez[match(names(maxgene.ct.ind),
                                                annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]}
    } else {
      print(paste0(rxn," NOT encoded by ct dependent gene!"))
    }   
}
##Show most abundant isoenzyme per ct
filter <- !apply(rxn.maxgprxct==0,1,any)
rxn.maxgprxct.f          <- rxn.maxgprxct[filter,]
rxn.maxgprxct.cellnote.f <- rxn.maxgprxct.cellnote[filter,]
rxn.repr.clusters.f      <- rxn.repr.clusters[filter,]
hm     <- heatmap.2(rxn.maxgprxct.f,symm = F,scale="none",notecol = "black",
                    cellnote = rxn.maxgprxct.cellnote.f,notecex=0.7,
                    col=brewer.pal(5,"Pastel2"),density.info="none",trace="none",
                    labRow=rxns.core.iso[match(rxn.repr.clusters.f,rxns.core.iso$RxnID),"Eq"])

#7: Save
writeLines(as.character(rxns.core.iso$RxnID),"rxns_core_iso.txt")
writeLines(as.character(rxns.core.pan$RxnID),"rxns_core_pan.txt")
writeLines(as.character(rxns.core.conn$RxnID),"rxns_core_conn.txt")

#Condense info
# #Describe subsystem coverage core rxns
# rxns.subs.min <- 50
# rxns.subs           <- table(HMR2.rxnTable[,"Subs"])
# rxns.core.subs      <- table(rxns.core[,"Subs"])/rxns.subs
# rxns.core.conn.subs <- table(rxns.core.conn[,"Subs"])/rxns.subs
# rxns.core.univ.subs <- table(rxns.core.univ[,"Subs"])/rxns.subs
# rxns.core.cont.subs <- table(rxns.core.cont[,"Subs"])/rxns.subs
# rxns.pan.fac.subs   <- table(rxns.pan.fac[,"Subs"])/rxns.subs
# rxns.pan.abs.subs   <- table(rxns.pan.abs[,"Subs"])/rxns.subs
# rxns.all.subs <-data.frame(rxns.core.conn.subs,
#                            rxns.core.cont.subs,
#                            rxns.core.univ.subs,
#                            rxns.pan.fac.subs,
#                            rxns.pan.abs.subs)
# rxns.all.subs.f              <- rxns.all.subs[rxns.subs>rxns.subs.min,]
# rxns.all.subs.sort           <- rxns.all.subs.f[order(rxns.core.subs[rxns.subs>rxns.subs.min],decreasing = T),]
# rownames(rxns.all.subs.sort) <- rxns.all.subs.sort[,1]
# rxns.all.subs.sort <- rxns.all.subs.sort[,-grep("Var1",colnames(rxns.all.subs.sort))]
# cols <- c("Red","Purple","Dark Red","Blue","Grey")
# barplot(as.matrix(t(rxns.all.subs.sort)),cex.names = 0.5,
#         beside = F,las=2,horiz=T,col=cols)