#Compare core-pan with housekeeping proteome, core-iso with tissue-enriched proteome,
#and cntx-clust with group enriched proteome
rm(list=ls())
setwd("~/Documents/Academia/Projects/RxnLands/Codes/analyzeGEMs/CompareToNormal")
require(VennDiagram)
require(gplots)
require(RColorBrewer)

#Load
load("../../data/annotation.RData")
load("../../data/preprocessed.RData") 
absn       <- read.delim("../RxnsCoreVsPan/rxns_absn.txt",stringsAsFactors=F)
  absn <- absn$RxnID
cntx_clust <- readLines("../RxnsCoreVsPan/rxns_cntx_clust.txt")
core_pan   <- readLines("../RxnsCoreVsPan/rxns_core_pan.txt")
core_iso   <- readLines("../RxnsCoreVsPan/rxns_core_iso.txt")
hpa.hk     <- read.delim("../../data/tissue_specificity_rna-any_expressed.tab",stringsAsFactors=F)
hpa.mix    <- read.delim("../../data/tissue_specificity_rna-any_Mixed.tab",stringsAsFactors=F)
hpa.all    <- read.delim("../../data/rna.csv",stringsAsFactors=T,sep = ",")
  HMR2.gpr   <- read.delim("../../docs/cHMR3765_rxnGeneMatrix.txt")
  rownames(HMR2.gpr)<-HMR2.gpr[,1]
HMR2.gpr    <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]
  geneIncMatFile    <- c("../geneInclMat.txt")
  geneInclMat       <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
HMR2.rxnTable <- read.delim("../../data/cHMR3765_rxnTable.txt",stringsAsFactors=F)
geneInclMat   <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
metgenes.EC   <- readLines("../../docs/ECdefinedMetGenes_HMR3765.txt")
  rxnIncMatFile     <- c("../rxnInclMat.txt")
  rxnInclMat        <- read.delim(rxnIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat  <- rxnInclMat[,-grep("X.1",colnames(rxnInclMat))]

#Param
percSampleWRxns          <- 0.95
arbitrary.threshold.rank <- 10
maxngpr.rxnint           <- 1
hpa.core.defin           <- "both"

#1: Preprocess
PSID.m        <- colnames(geneInclMat)
  lib.size          <- with(y$samples, lib.size * norm.factors)
  E                 <- t(log2(t(y$counts + 0.5)/(lib.size + 1) * 1e+06))
  E.ensg            <- E
  colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
  rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat)) #This is the reference metabolic genome
geneInclMat.c     <- geneInclMat[genes.ensg.common,]
genes.ensg.m      <- rownames(geneInclMat.c[which(apply(geneInclMat.c!=0,1,any)),])
  E.ensg.c          <- E.ensg[genes.ensg.common,] #n common genes x nS (1082)
  PSID         <- rownames(y$samples)
  PSID.s       <- gsub(pattern="-",replacement=".",x=PSID)
  filter       <- PSID.s%in%colnames(geneInclMat)
  mapper       <- match(colnames(geneInclMat),PSID.s[filter])
E.ensg.c.m   <- E.ensg.c[,filter][,mapper] #n common genes x nS in models (917)

###Above: the reference genome for the reconstructed models

#3: Classify HPA proteins
hpa.tissues    <- levels(hpa.all$Sample)[!grepl("[A-Z]",
                                                levels(hpa.all$Sample))]
hpa.onlytissue <- hpa.all[hpa.all$Sample %in% hpa.tissues,]
hpa.abundxgene <- aggregate(Abundance~Gene,data=hpa.onlytissue,FUN = table,simplify = T)
min.ntissues   <- length(hpa.tissues)*(1-percSampleWRxns)
hpa.core.genes  <- hpa.abundxgene$Gene[hpa.abundxgene$Abundance
                                      [,"Not detected"]<min.ntissues]
  hpa.core.genes <- as.character(hpa.core.genes)
hpa.cntx.genes <- hpa.abundxgene$Gene[hpa.abundxgene$Abundance
                                      [,"Not detected"]>min.ntissues &
                                       hpa.abundxgene$Abundance
                                      [,"Not detected"]<length(hpa.tissues)]
  hpa.cntx.genes <- as.character(hpa.cntx.genes)
hpa.absn.genes <- hpa.abundxgene$Gene[hpa.abundxgene$Abundance
                                      [,"Not detected"]==length(hpa.tissues)]
  hpa.absn.genes <- as.character(hpa.absn.genes)

###Go metabolic
if (hpa.core.defin != "my"){#Use HPA or my or both way to define housekeeping?
  hpa.hk.genes   <- hpa.hk$Ensembl[hpa.hk$Ensembl %in% genes.ensg.common]
  hpa.mix.genes  <- hpa.mix$Ensembl[hpa.mix$Ensembl %in% genes.ensg.common]
  hpa.core.genes.hpa  <- unique(c(hpa.hk.genes,hpa.mix.genes))
  if (hpa.core.defin == "both"){
    hpa.core.genes.my  <- hpa.core.genes[hpa.core.genes %in% genes.ensg.common]
    hpa.core.genes     <- unique(c(hpa.core.genes.hpa,hpa.core.genes.my))
  } else {hpa.core.genes <- hpa.core.genes.hpa}} else {
    hpa.core.genes  <- hpa.core.genes[hpa.core.genes %in% genes.ensg.common]
  }
hpa.absn.genes <- as.character(hpa.absn.genes[hpa.absn.genes %in% genes.ensg.common])
hpa.cntx.genes <- as.character(hpa.cntx.genes[hpa.cntx.genes %in% genes.ensg.common])
  hpa.cntx.genes <- hpa.cntx.genes[!hpa.cntx.genes %in% hpa.core.genes]
hpa.unch.genes <- setdiff(genes.ensg.common,c(hpa.absn.genes,hpa.core.genes,hpa.cntx.genes))
  
#4: Classify reactions according to HPA
hpa.core.rxns <- rownames(HMR2.gpr)[apply(HMR2.gpr[,hpa.core.genes]!=0,1,any)]
hpa.cntx.rxns <- rownames(HMR2.gpr)[apply(HMR2.gpr[,hpa.cntx.genes]!=0,1,any)]  
hpa.absn.rxns <- rownames(HMR2.gpr)[apply(HMR2.gpr[,hpa.absn.genes]!=0,1,any)]  
hpa.unch.rxns <- rownames(HMR2.gpr)[apply(HMR2.gpr[,hpa.unch.genes]!=0,1,any)]
ncheckedrxns  <- sum(length(hpa.core.rxns)+length(hpa.cntx.rxns)+length(hpa.absn.rxns))

writeLines(hpa.core.rxns,"hpa_rxns_core.txt",sep = "\n")
writeLines(hpa.cntx.rxns,"hpa_rxns_cntx.txt",sep = "\n")
writeLines(hpa.absn.rxns,"hpa_rxns_absn.txt",sep = "\n")
writeLines(hpa.unch.rxns,"hpa_rxns_unch.txt",sep = "\n")

#5: Check how many core_pan rxns are core in HPA
core <- c(core_pan,core_iso) #Core is core, regardless of gene assoc
core <- core[-which(core=="HMR_8426")] #ARG1 is not core - just misannotated
only.core          <- sum(!core %in% hpa.core.rxns)
only.hpa.core      <- sum(!hpa.core.rxns %in% core)
both.core.hpa.core <- sum(core %in% hpa.core.rxns)
not.core.hpa.core  <- ncheckedrxns - both.core.hpa.core
                                  - only.core - only.hpa.core
grid.newpage()
venn.plot   <- draw.pairwise.venn(length(core),length(hpa.core.rxns),
                                  both.core.hpa.core,category=c("Core - pan",
                                                                    "Housekeeping\nproteome"),fill=brewer.pal(11,"RdYlBu")[c(1,11)],
                                  cat.pos=c(315,45),cat.dist=0.06,cat.cex=2,
                                  cex=2)
grid.draw(venn.plot)

cont.table  <- matrix(c(both.core.hpa.core,only.core,
                        only.hpa.core,not.core.hpa.core),ncol=2,nrow=2)
fisher.test <- fisher.test(cont.table)

# #4: Find genes of core_pan
# core_pan.genes <- colnames(HMR2.gpr)[apply(HMR2.gpr[core_pan,]!=0,2,any)]
#   core_pan.genes <- core_pan.genes[core_pan.genes %in% genes.ensg.m]
# 
# only.core_pan          <- sum(!core_pan.genes %in% hpa.core.genes)
# only.hpa.core          <- sum(!hpa.core.genes %in% core_pan.genes)
# both.core_pan.hpa.core <- sum(core_pan.genes %in% hpa.core.genes)
# not.core_pan.hpa.core  <- length(genes.ensg.common) - both.core_pan.hpa.core
#                         - only.core_pan - only.hpa.core
# grid.newpage()
# venn.plot   <- draw.pairwise.venn(length(core_pan.genes),length(hpa.core.genes),
#                                 both.core_pan.hpa.core,category=c("Core - pan",
#                                 "Housekeeping\nproteome"),fill=brewer.pal(11,"RdYlBu")[c(1,11)],
#                                 cat.pos=c(315,45),cat.dist=0.06,cat.cex=2,
#                                 cex=2)
# grid.draw(venn.plot)
# 
# cont.table  <- matrix(c(both.core_pan.hpa.core,only.core_pan,
#                        only.hpa.core,not.core_pan.hpa.core),ncol=2,nrow=2)
# fisher.test <- fisher.test(cont.table)

###Which rxns are only core and not hk?
only.core.rxns <- core[!core %in% hpa.core.rxns]
  becauseUnchecked <- only.core.rxns %in% hpa.unch.rxns
  becauseCntx      <- only.core.rxns %in% hpa.cntx.rxns
  becauseAbsn      <- only.core.rxns %in% hpa.absn.rxns #If none of these, then the genes are not in the "common" stack
only.core.rxns.tocheck <- only.core.rxns[becauseAbsn |
                                         becauseCntx |
                                         becauseUnchecked]
writeLines(only.core.rxns.tocheck,"rxns_core_onlyincancer.txt")

# ###Select most interesting: 1to1 gpr, highly expressed, low variability, EC
#   only.core_pan.genes.ngpr   <- apply(HMR2.gpr[,only.core_pan.genes]!=0,2,sum)
#   only.core_pan.genes.1gpr   <- names(only.core_pan.genes.ngpr)[only.core_pan.genes.ngpr<=maxngpr.rxnint]
#   only.core_pan.genes.1gpr.EC   <- only.core_pan.genes.1gpr[only.core_pan.genes.1gpr %in% metgenes.EC]
#   E.only.core_pan.genes.1gpr.EC <- E.ensg.c.m[rownames(E.ensg.c.m) %in% only.core_pan.genes.1gpr.EC,]
#   var.score <- rank(apply(E.only.core_pan.genes.1gpr.EC,1,mad))
#   exp.score <- rank(-apply(E.only.core_pan.genes.1gpr.EC,1,mean))
#   tot.score <- 0.5*var.score + 0.5*exp.score
# only.core_pan.genes.top.score <- tot.score[order(tot.score)]
# only.core_pan.genes.top       <- names(only.core_pan.genes.top.score)[1:arbitrary.threshold.rank]

###How are their associated gene expressed in the models?
only.core.genes <- c()
for (rxn in only.core.rxns.tocheck){
  rxn.ex.genes       <- colnames(HMR2.gpr)[apply(HMR2.gpr[rxn,]!=0,2,any)]
    rxn.ex.genes     <- rxn.ex.genes[rxn.ex.genes %in% genes.ensg.common]
  rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
  E.ensg.c.zerneg    <- E.ensg.c.m
    E.ensg.c.zerneg[E.ensg.c.m<0] <- 0
  rxn.ex.gene.sort   <- sort.int(E.ensg.c.zerneg[rxn.ex.genes[1],],
                                          index.return = T)
  cols5              <- brewer.pal(5,"Set2")
  barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
          col=cols5,legend.text = ifelse(rxn.ex.genes.hgnc=="",rxn.ex.genes,rxn.ex.genes.hgnc),
         horiz=T,
          main=HMR2.rxnTable[HMR2.rxnTable[,"RxnID"]==rxn,"Eq"],
          xlab="Size-corrected log-cpm")
  only.core.genes <- c(only.core.genes,rxn.ex.genes)}
only.core.genes <- unique(only.core.genes)
only.core.genes <- only.core.genes[-which(only.core.genes=="ENSG00000118520")]
writeLines(only.core.genes,"genes_core_onlyincancer.txt")

#5bis: Show heatmap of core and hpa detection
rnxs.coreorabs <- only.core.rxns.tocheck
rxns.coreorabs.eq <- HMR2.rxnTable[match(rnxs.coreorabs,HMR2.rxnTable$RxnID),"Eq"] 
gene.coreorabs <- only.core.genes
hmcols  <- colorRampPalette(brewer.pal(9, "Blues"))(256)
hm      <- heatmap.2(E.ensg.c.m[gene.coreorabs,],Rowv = T,Colv = T,symm = F,
                     scale="none",col=hmcols,density.info="none",trace="none",
                     labCol=NA,labRow = annotEntrez[match(gene.coreorabs,annotEntrez[,"ensembl_gene_id"])
                                                    ,"hgnc_symbol"],dendrogram="col",breaks=seq(-1,10,by=11/256))
geneInclMat.c.ord <- as.matrix(geneInclMat.c[gene.coreorabs,]
                               [rev(hm$rowInd),hm$colInd])
hm.2    <- heatmap.2(geneInclMat.c.ord,Rowv = F,Colv = F,symm = F,
                     scale="none",col=c("White","Black"),density.info="none",trace="none",
                     labCol=NA,labRow = annotEntrez[match(gene.coreorabs[rev(hm$rowInd)],annotEntrez[,"ensembl_gene_id"])
                                                    ,"hgnc_symbol"],dendrogram="none")
rxnInclMat.ord    <- as.matrix(rxnInclMat[rxns.coreorabs.eq,hm$colInd])
hm.3    <- heatmap.2(rxnInclMat.ord,Rowv = F,Colv = F,symm = F,
                     scale="none",col=c("White","Black"),density.info="none",trace="none",
                     labCol=NA,labRow = rxns.coreorabs.eq,
                     dendrogram="none")
hpa.abund.gene.coreorabs      <- hpa.abundxgene$Abundance[match(gene.coreorabs,hpa.abundxgene$Gene),]
rownames(hpa.abund.gene.coreorabs) <- annotEntrez[match(gene.coreorabs,
                                                        annotEntrez[,"ensembl_gene_id"]),
                                                  "hgnc_symbol"]
hpa.abund.gene.coreorabs.ord  <- hpa.abund.gene.coreorabs[rev(hm$rowInd),
                                                          c("Not detected","Low","Medium","High")]
hpa.abund.gene.coreorabs.frac <- hpa.abund.gene.coreorabs.ord/length(hpa.tissues)
barplot(t(hpa.abund.gene.coreorabs.frac),legend.text = T,las=2,
        col=c("white",brewer.pal(3,"YlOrBr")),ylab="Fraction of human tissues")

#6: Check how many absent rxns are core in HPA
only.absn          <- sum(!absn %in% hpa.core.rxns)
only.hpa.core      <- sum(!hpa.core.rxns %in% absn)
both.absn.hpa.core <- sum(absn %in% hpa.core.rxns)
not.absn.hpa.core  <- ncheckedrxns - both.absn.hpa.core
                      - only.absn - only.hpa.core
grid.newpage()
venn.plot   <- draw.pairwise.venn(length(absn),length(hpa.core.rxns),
                                  both.absn.hpa.core,category=c("Absent",
                                  "Housekeeping\nproteome"),
                                  fill=brewer.pal(11,"RdYlBu")[c(7,11)],
                                  cat.pos=c(315,45),cat.dist=0.06,cat.cex=2,
                                  cex=2)
grid.draw(venn.plot)

cont.table  <- matrix(c(both.absn.hpa.core,only.absn,
                        only.hpa.core,not.absn.hpa.core),ncol=2,nrow=2)
fisher.test <- fisher.test(cont.table)

###Which rxns are only absent and not hk?
only.absn.rxns     <- absn[absn %in% hpa.core.rxns]

# #4: Find genes of absn
# absn.genes  <- rownames(geneInclMat.c[which(apply(geneInclMat.c==0,1,all)),])
# 
# only.absn          <- sum(!absn.genes %in% hpa.core.genes)
# only.hpa.core      <- sum(!hpa.core.genes %in% absn.genes)
# both.absn.hpa.core <- sum(absn.genes %in% hpa.core.genes)
# not.absn.hpa.core  <- length(genes.ensg.common) - both.absn.hpa.core
#                                 - only.absn - only.hpa.core
# grid.newpage()
# venn.plot   <- draw.pairwise.venn(length(absn.genes),length(hpa.core.genes),
#                                   both.absn.hpa.core,category=c("Absent",
#                                     "Housekeeping\nproteome"),
#                                   fill=brewer.pal(11,"RdYlBu")[c(7,11)],
#                                   cat.pos=c(315,45),cat.dist=0.06,cat.cex=2,
#                                   cex=2)
# grid.draw(venn.plot)
# 
# cont.table  <- matrix(c(both.absn.hpa.core,only.absn,
#                         only.hpa.core,not.absn.hpa.core),ncol=2,nrow=2)
# fisher.test <- fisher.test(cont.table)
# 
# ###Which rxns are only absn but hk?
# only.absn.genes <- absn.genes[absn.genes %in% hpa.core.genes]
# writeLines(only.absn.genes,"absn_genesbuthpapan.txt")

###Select most interesting: 1to1 gpr, lowly expressed, low variability, EC, high expressed in hpa
only.absn.rxns.ngpr       <- apply(HMR2.gpr[only.absn.rxns,]!=0,1,sum)
only.absn.rxns.1gpr       <- names(only.absn.rxns.ngpr)[only.absn.rxns.ngpr<=maxngpr.rxnint]
only.absn.genes.1gpr      <- colnames(HMR2.gpr)[apply(HMR2.gpr[only.absn.rxns.1gpr,]!=0,2,any)]
only.absn.genes.1gpr.EC   <- only.absn.genes.1gpr[only.absn.genes.1gpr %in% metgenes.EC]
E.only.absn.genes.1gpr.EC <- E.ensg.c.m[rownames(E.ensg.c.m) %in% only.absn.genes.1gpr.EC,]
  var.score <- rank(apply(E.only.absn.genes.1gpr.EC,1,mad))
  exp.score <- rank(apply(E.only.absn.genes.1gpr.EC,1,mean))
  hpa.scores.all <- apply(sweep(hpa.abundxgene$Abundance,2,c(10,1,5,0),"*"),1,sum)
  hpa.scores.genes.1gpr.EC  <- hpa.scores.all[match(only.absn.genes.1gpr.EC,hpa.abundxgene$Gene)]
  hpa.scores     <- rank(-hpa.scores.genes.1gpr.EC)
tot.score <- 0.05*var.score + 0.15*exp.score + 0.8*hpa.scores

only.absn.genes.top.score <- tot.score[order(tot.score)]
only.absn.genes.top       <- names(only.absn.genes.top.score)[1:arbitrary.threshold.rank]

###Plot top
only.absn.genes.top.rxns <- c()
for (gene.top in only.absn.genes.top){
  only.absn.rxns   <- rownames(HMR2.gpr)[which(HMR2.gpr[,gene.top]!=0)]
  only.absn.rxns   <- only.absn.rxns[only.absn.rxns %in% absn]
  E.ensg.c.zerneg               <- E.ensg.c.m
  E.ensg.c.zerneg[E.ensg.c.m<0] <- 0
  rxn.ex.genes.hgnc           <- annotEntrez[match(gene.top,
                                                   annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
  rxn.ex.gene.sort            <- sort.int(E.ensg.c.zerneg[gene.top,],
                                          index.return = T)
  cols5                       <- brewer.pal(5,"Set2")
  barplot(E.ensg.c.zerneg[gene.top,rxn.ex.gene.sort$ix],beside=F,las=2,
          col=cols5,legend.text = ifelse(rxn.ex.genes.hgnc=="",gene.top,rxn.ex.genes.hgnc),
          horiz=T,
          main=HMR2.rxnTable[HMR2.rxnTable[,"RxnID"]==only.absn.rxns,"Eq"],
          xlab="Size-corrected log-cpm",names.arg=NA)
  only.absn.genes.top.rxns <- c(only.absn.genes.top.rxns,only.absn.rxns)}
only.absn.genes.top.rxns <- unique(only.absn.genes.top.rxns)

#7: Show heatmap of top absent - pan and hpa detection
rnxs.coreorabs <- c(only.core.rxns.tocheck,only.absn.genes.top.rxns)
rxns.coreorabs.eq <- HMR2.rxnTable[match(rnxs.coreorabs,HMR2.rxnTable$RxnID),"Eq"] 
gene.coreorabs <- c(only.core.genes,only.absn.genes.top)
hmcols  <- colorRampPalette(brewer.pal(9, "Blues"))(256)
hm      <- heatmap.2(E.ensg.c.m[gene.coreorabs,],Rowv = T,Colv = T,symm = F,
                     scale="none",col=hmcols,density.info="none",trace="none",
                     labCol=NA,labRow = annotEntrez[match(gene.coreorabs,annotEntrez[,"ensembl_gene_id"])
                         ,"hgnc_symbol"],dendrogram="col",breaks=seq(-1,10,by=11/256))
geneInclMat.c.ord <- as.matrix(geneInclMat.c[gene.coreorabs,]
                               [hm$rowInd,hm$colInd])
hm.2    <- heatmap.2(geneInclMat.c.ord,Rowv = F,Colv = F,symm = F,
                     scale="none",col=c("White","Black"),density.info="none",trace="none",
                     labCol=NA,labRow = annotEntrez[match(gene.coreorabs[hm$rowInd],annotEntrez[,"ensembl_gene_id"])
                       ,"hgnc_symbol"],dendrogram="none")
rxnInclMat.ord    <- as.matrix(rxnInclMat[rxns.coreorabs.eq,hm$colInd])
hm.3    <- heatmap.2(rxnInclMat.ord,Rowv = F,Colv = F,symm = F,
                     scale="none",col=c("White","Black"),density.info="none",trace="none",
                     labCol=NA,labRow = rxns.coreorabs.eq,
                     dendrogram="none")
hpa.abund.gene.coreorabs      <- hpa.abundxgene$Abundance[match(gene.coreorabs,hpa.abundxgene$Gene),]
  rownames(hpa.abund.gene.coreorabs) <- annotEntrez[match(gene.coreorabs,
                                        annotEntrez[,"ensembl_gene_id"]),
                                        "hgnc_symbol"]
hpa.abund.gene.coreorabs.ord  <- hpa.abund.gene.coreorabs[hm$rowInd,
                                    c("Not detected","Low","Medium","High")]
hpa.abund.gene.coreorabs.frac <- hpa.abund.gene.coreorabs.ord/length(hpa.tissues)
barplot(t(hpa.abund.gene.coreorabs.frac),legend.text = T,las=2,
        col=c("white",brewer.pal(3,"YlOrBr")),ylab="Fraction of human tissues")

#8: Check levels of the cntx rxns in HPA
tissues        <- c("cerebral cortex","cerebral cortex","kidney","pancreas",
                    "rectum","colon","endometrium","ovary","lung",
                    "urinary bladder","esophagus","breast","lung")
cntx_clust.hpamaxexpr <- matrix(0,nrow=length(cntx_clust),ncol=length(tissues),
                                dimnames=list(cntx_clust,tissues))
for (rxn in cntx_clust){
  rxn.genes     <- colnames(HMR2.gpr)[apply(HMR2.gpr[rxn,]!=0,2,any)]
  rxn.genes     <- rxn.genes[rxn.genes %in% genes.ensg.common]
  hpa.genes     <- hpa.onlytissue[hpa.onlytissue$Gene %in% rxn.genes,]
  if (nrow(hpa.genes)!=0){
    hpa.maxfpkm   <- aggregate(Value~Sample,data=hpa.genes,FUN = max,simplify = T)
    cntx_clust.hpamaxexpr[rxn,] <- hpa.maxfpkm[match(tissues,hpa.maxfpkm$Sample),"Value"]
  } else {cntx_clust.hpamaxexpr[rxn,] <- NA}
}
rownames(cntx_clust.hpamaxexpr) <- HMR2.rxnTable[match(cntx_clust,
                                                 HMR2.rxnTable$RxnID),"Eq"]
colnames(cntx_clust.hpamaxexpr) <- toupper(tissues)
hmcols  <- colorRampPalette(c("white",brewer.pal(9, "YlOrRd")))(256)
hm      <- heatmap.2(cntx_clust.hpamaxexpr,Rowv = F,Colv = F,symm = F,
                     scale="none",density.info="none",trace="none",
                     col=hmcols,breaks=seq(0,50,by=50/256))

#9: Define cntx reactions unique to tumors
cntx_clust.uniqtoT.eq <- c("mannose[s] + Na+[s] => mannose[c] + Na+[c]",
                           "H2O[c] + N-acetylneuraminate-9-phosphate[c] => N-acetylneuraminate[c] + Pi[c]",
                           "cholesterol[c] <=> cholesterol[s]",
                           "IV2Fuc-Lc4Cer[c] + UDP-N-acetyl-D-galactosamine[c] => type I A glycolipid[c] + UDP[c]",                                                                                
                           "IV2Fuc-Lc4Cer[c] + UDP-galactose[c] => type I B glycolipid[c] + UDP[c]")
cntx_clust.uniqtoT    <- HMR2.rxnTable[HMR2.rxnTable$Eq %in% cntx_clust.uniqtoT.eq,
                                       "RxnID"]
only.cntx.genes <- c()
for (rxn in cntx_clust.uniqtoT){
  rxn.genes     <- colnames(HMR2.gpr)[apply(HMR2.gpr[rxn,]!=0,2,any)]
  rxn.genes     <- rxn.genes[rxn.genes %in% genes.ensg.common]
  only.cntx.genes <- c(only.cntx.genes,rxn.genes)
}
writeLines(only.cntx.genes,"../../correlatePARIS/data/only.cntx.genes.txt")
