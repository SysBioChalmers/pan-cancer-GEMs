#Enrich mets in class rxns against HMR2
rm(list=ls())
require(synapseClient)
synapseLogin("gatto@chalmers.se","forre$ter1867")

#Param
significance.level <- 0.01
class              <- "core_iso"

#Load
rxntable.file <- synGet("syn3218826")
rxnTable      <- read.delim(rxntable.file@filePath,row.names=1,stringsAsFactors=F)
topSub        <- read.delim(paste0("rxns_",class,".txt"),sep="\t",header=F)
nmetsinnet    <- read.delim(paste0("nrxnmetin_rxns_",class,"vsHMR2.txt"),sep="\t",header=T)

#1: Get counts
topSub.rxns    <- topSub[,"V1"]
topSub.nrxns   <- length(topSub.rxns)
HMR3765.nrxns  <- dim(rxnTable)[1]
HMR3765.subs   <- nmetsinnet[,"X"]
rownames(nmetsinnet)  <- HMR3765.subs

#2: Fisher test
sub.enrichment <- matrix(0,nrow=length(HMR3765.subs),ncol=2,
                         dimnames=list(HMR3765.subs,c("p-value","enrichment")))
for (sub in HMR3765.subs){
  topSub.nrxns.sub    <- nmetsinnet[sub,grep("InList",colnames(nmetsinnet))]
  topSub.nrxns.notsub <- topSub.nrxns - topSub.nrxns.sub
  HMR3765.nrxns.sub   <- nmetsinnet[sub,grep("InHMR2",colnames(nmetsinnet))]
  HMR3765.nrxns.nosub <- (HMR3765.nrxns-topSub.nrxns) - HMR3765.nrxns.sub
  
  contTable = matrix(c(topSub.nrxns.sub,topSub.nrxns.notsub,
                       HMR3765.nrxns.sub,HMR3765.nrxns.nosub),nrow=2,ncol=2)
  fisher <- fisher.test(contTable)
  sub.enrichment[sub,] <- c(fisher$p.value,fisher$estimate)}

#3: Sort results: filter enrichment that is below significance level (after adj for multiple testing)
j                     <- sort.int(sub.enrichment[,"enrichment"],decreasing=F,index.return=T)
adj.pvalue            <- p.adjust(sub.enrichment[,"p-value"],method="fdr")
sub.enrichment.annot  <- cbind(sub.enrichment,nmetsinnet[,c(grep("InHMR2",colnames(nmetsinnet)),grep("InList",colnames(nmetsinnet)))])
sub.enrichment.sort   <- sub.enrichment.annot[j$ix,]
adj.pvalue            <- p.adjust(sub.enrichment.sort[,"p-value"],method="fdr")
sub.enrichment.adj    <- cbind(sub.enrichment.sort,adj.pvalue)
sub.enrichment.sign   <- sub.enrichment.adj[sub.enrichment.adj[,"adj.pvalue"]<significance.level,]
sub.enrichment.sign.under <- sub.enrichment.sign[sub.enrichment.sign[,"enrichment"]<1,]
sub.enrichment.sign.over  <- sub.enrichment.sign[sub.enrichment.sign[,"enrichment"]>1,]

#4: Plot
barplot(log2(sub.enrichment.sign[,"enrichment"]+2^-6),ylim=c(0,6),
        names.arg=rownames(sub.enrichment.sign),
        ylab=paste("Enrichment of metabolites in",class,"rxns"),
        col=c(rep("dark blue",dim(sub.enrichment.sign.under)[1]),
                  rep("dark red",dim(sub.enrichment.sign.over)[1])),las=2)