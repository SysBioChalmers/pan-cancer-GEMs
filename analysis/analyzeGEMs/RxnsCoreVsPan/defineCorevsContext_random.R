# Pre stuff - adapted from HMR_stats.R by L. Väremo (2015):
# setwd("~/Box Sync/3rd Semester (May 2013 - November 2013)/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan")
require(gplots)
library(openxlsx)
require(RColorBrewer)

# Parameters
fracModels  <- 0.95
nRandModels <- 917
nPerms      <- 1:1000

# Read HMR from Excel file:
hmr <- read.xlsx("HMR3765_20140415.xlsx",sheet="RXNS",detectDates=F)
hmr <- hmr[is.na(hmr[,1]),] # Remove commented (#) lines (non-commented should be <NA>s)
hmr <- hmr[,colnames(hmr)%in%c("RXNID","EQUATION","GENE.ASSOCIATION","SUBSYSTEM")]

# Get nGenes distribution in real models
geneIncMatFile <- c("geneInclMat.txt")
geneInclMat    <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
geneInclMat    <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
nGenes.xmodel  <- apply(geneInclMat==1,2,sum)

# Core/pan reactions for different number of genes:
genepool <- unique(unlist(strsplit(hmr$GENE.ASSOCIATION,";")))
rxn2gene <- strsplit(hmr$GENE.ASSOCIATION,";")
myFun    <- function(nPerms,genepool,rxn2gene,nGenes.xmodel,nRandModels,fracModels) {
             core <- rep(NA,length(nPerms))
             absn <- rep(NA,length(nPerms))
             cntx <- rep(NA,length(nPerms))
             size <- matrix(ncol=4,nrow=length(nPerms))
             for(j in 1:length(nPerms)) {
                nGenes.xmodel.rand <- sample(nGenes.xmodel,nRandModels)
                rxnmat <- matrix(ncol=nRandModels,nrow=length(rxn2gene))
                for(i in 1:nRandModels) {
                  print(paste("Doing",i))
                   genesample <- sample(genepool,nGenes.xmodel.rand[i])
                   rxnmat[,i] <- unlist(lapply(rxn2gene,function(x) ifelse(!all(!x%in%genesample),1,0)))
                }
                core[j]  <- sum(apply(rxnmat,1,sum)>=nRandModels*fracModels) # reactions in X% of "models"
                absn[j]  <- sum(apply(rxnmat,1,sum)<=nRandModels*(1-fracModels)) # reactions in (1-X)% of "models"
                cntx[j]  <- sum(apply(rxnmat,1,sum)>nRandModels*(1-fracModels) & apply(rxnmat,1,sum)<nRandModels*fracModels) # reactions between (1-X)% and X% of "models"                
                size[j,] <- summary(apply(rxnmat,2,sum))[c(1,3,4,6)]
             }
             stats <- cbind(nPerms,core,absn,cntx,size)
}
#
# Run on Glenn for faster runtime...
library(snowfall)
sfInit(parallel=T,cpus=16)
out <- sfLapply(nPerms,myFun,genepool,rxn2gene,nGenes.xmodel,nRandModels,fracModels) # this will be a list, need to remake to a matrix...
stats <- matrix(unlist(out),ncol=8,byrow=T)
#...or run locally:
# stats <- myFun(nPerms,genepool,rxn2gene)
#
colnames(stats) <- c("nPerms","nCoreRxns","nAbsnRxns","nCntxRxns","minRxns","medianRxns","meanRxns","maxRxns")
#
# barplot(rbind(stats[,2],stats[,5]-stats[,2]),names.arg=size[,1],las=2,ylim=c(0,8000),
#         main="Mean and core reactions",ylab="Number of reactions",xlab="Number of genes")
save("stats",file="ncorevscntx_randommodels.RData")

#5: Compute stats on random core vs. absn vs. cntx
load("ncorevscntx_randommodels.RData")
quantiles <- apply(stats,2,quantile,seq(0,1,1/20))
n   <- quantiles["50%",c("nCoreRxns","nCntxRxns","nAbsnRxns")] 
n.l <- quantiles["5%",c("nCoreRxns","nCntxRxns","nAbsnRxns")] 
n.u <- quantiles["95%",c("nCoreRxns","nCntxRxns","nAbsnRxns")] 
barplot2(n,horiz=F,beside=T,las=2,names.arg = c("Core","Contextual","Absent"),         
         ci.l =n.l,ci.u =n.u,
         col = brewer.pal(length(n),"Greys")[c(3,2,1)],
         ylab="N° reactions",plot.ci=T,
         legend.text=F,ylim=c(0,7000))