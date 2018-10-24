###Define core vs. contextual rxns in cancer GEMs vs. normal GEMs
rm(list=ls())
setwd("~/Documents/Academia/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan")
require(boot)
library(gplots)
library(ggplot2)
require(RColorBrewer)
require(vioplot)

#Load: gene expression, models (in terms of rxn and gene inclusion matrices), HMR2 gpr and rxns
load("../../data/preprocessed.RData") 
  rxnIncMatFile <- c("../rxnInclMat.txt")
  rxnInclMat    <- read.delim(rxnIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
rxnInclMat    <- rxnInclMat[,-grep("X.1",colnames(rxnInclMat))]
  geneIncMatFile <- c("../geneInclMat.txt")
  geneInclMat    <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
geneInclMat    <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
HMR2.rxnTable <- read.delim("../../docs/cHMR3765_rxnTable.txt")
  HMR2.gpr  <- read.delim("../../docs/cHMR3765_rxnGeneMatrix.txt")
  rownames(HMR2.gpr)<-HMR2.gpr[,1]
HMR2.gpr <- HMR2.gpr[,-c(1,ncol(HMR2.gpr))]

#Param
nCancerTypes    <- 13
percSampleWRxns <- 0.95
cols            <- colorRampPalette(brewer.pal(12, "Set3"))(nCancerTypes)

#1: Preprocess
### Assign each model a cancer type
PSID         <- rownames(y$samples)
PSID.s       <- gsub(pattern="-",replacement=".",x=PSID)
filter       <- PSID.s%in%colnames(rxnInclMat)
mapper       <- match(colnames(rxnInclMat),PSID.s[filter])
fGlobalCT.m  <- fGlobalCT[filter][mapper]
names(fGlobalCT.m)  <- PSID.s[filter][mapper]

#2: How many n° rxns per model?
GEMS.rxnXPSID <- apply(rxnInclMat==1,2,sum) #rxnInclMat=1 indicates presence of the rxn in the model, 0 otherwise
par(mfrow=c(1,2))
boxplot2(GEMS.rxnXPSID ~ fGlobalCT.m,col=cols,ylim=c(0,max(GEMS.rxnXPSID)),
        ylab="N° reactions in the model",las=2)
vioplot(GEMS.rxnXPSID,col="Grey",names = "ALL",ylim=c(0,max(GEMS.rxnXPSID)))
par(mfrow=c(1,1))

#2bis: How many n° genes per model?
GEMS.geneXPSID <- apply(geneInclMat==1,2,sum) #geneInclMat=1 indicates presence of the rxn in the model, 0 otherwise
par(mfrow=c(1,2))
boxplot2(GEMS.geneXPSID ~ fGlobalCT.m,col=cols,ylim=c(0,max(GEMS.geneXPSID)),
         ylab="N° genes in the model",las=2,names=toupper(levels(fGlobalCT.m)))
vioplot(GEMS.geneXPSID,col="Grey",names = "ALL",ylim=c(0,max(GEMS.geneXPSID)))
par(mfrow=c(1,1))

#3: How many models have a given rxn?
nModelsXRxn      <- apply(rxnInclMat==1,1,sum)
nModelsXRxn.perc <- nModelsXRxn/ncol(rxnInclMat)
barplot(nModelsXRxn.perc[order(nModelsXRxn.perc)],names.arg = NA,border=NULL,
        xlab="HMR2 reactions", ylab="Fraction of models with the reaction")

# #3bis: When does the fraction of models with a given rxn flattens?
# l          <- lowess(nModelsXRxn.perc[order(nModelsXRxn.perc)],f=1/50)
# lines(l,col="green")
# derivative <- l$x
# slope.diff <- l$x
# for (i in 2:length(nModelsXRxn.perc)){
#   derivative[i] <- (l$y[i]-l$y[i-1])/(l$x[i]-l$x[i-1])
#   if (i==2){next}
#   slope.diff[i] <- abs(derivative[i]-derivative[i-1])}
# slope.diff[1:length(nModelsXRxn.perc)*0.9] <-0
# nModelsXRxn.perc.opt <- l$y[which.max(slope.diff)]
# abline(h=nModelsXRxn.perc.opt,col="red")

#4: How many rxns are declared core or context given a threshold for presence in x% of the models?
percsSampleWRxns <- c(seq(0.5,0.90,by = 0.05),
                      seq(0.91,0.99,by = 0.01),
                      seq(0.991,1,by = 0.001))
nRxnsClassified  <- matrix(0,nrow=3,ncol=length(percsSampleWRxns),
                           dimnames=list(c("N°core","N°context","N°absent"),
                                         paste0("Core@",percsSampleWRxns)))
for (perc in percsSampleWRxns){
  rxns.core       <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc>perc]),])
  rxns.context    <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc<=perc&nModelsXRxn.perc>0]),])
  rxns.abs        <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc==0]),])
  nRxnsClassified[,paste0("Core@",perc)] <- rbind(nrow(rxns.core),
                                                  nrow(rxns.context),
                                                  nrow(rxns.abs))}
barplot(nRxnsClassified,beside = F,las=2,ylab="N° rxns per group")

#5: Establish classification - Let's settle for 0.95
rxns.core <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc>percSampleWRxns]),],"Type"="Core")
rxns.cntx <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc<=percSampleWRxns&nModelsXRxn.perc>0]),],"Type"="Pan","Class"="Facultative")
rxns.absn <- data.frame(HMR2.rxnTable[HMR2.rxnTable[,"Eq"]%in%names(nModelsXRxn.perc[nModelsXRxn.perc==0]),],"Type"="Pan","Class"="Absent")

#6: Bootstrap statistics
nbootstraps <- 1000
n.core    <- function(modelsxrxn,i) sum(modelsxrxn[i]>percSampleWRxns)
n.core.bs <- boot(nModelsXRxn.perc,statistic = n.core,R = nbootstraps)
n.core.ci <- boot.ci(n.core.bs,type="basic")
n.cntx    <- function(modelsxrxn,i) sum(modelsxrxn[i]<percSampleWRxns & modelsxrxn[i]>0)
n.cntx.bs <- boot(nModelsXRxn.perc,statistic = n.cntx,R = nbootstraps)
n.cntx.ci <- boot.ci(n.cntx.bs,type="basic")
n.absn    <- function(modelsxrxn,i) sum(modelsxrxn[i]==0)
n.absn.bs <- boot(nModelsXRxn.perc,statistic = n.absn,R = nbootstraps)
n.absn.ci <- boot.ci(n.absn.bs,type="basic")

n   <- c(n.core.bs$t0,n.cntx.bs$t0,n.absn.bs$t0)
n.l <- c(n.core.ci$basic[1,4],n.cntx.ci$basic[1,4],n.absn.ci$basic[1,4])
n.u <- c(n.core.ci$basic[1,5],n.cntx.ci$basic[1,5],n.absn.ci$basic[1,5])

barplot2(n,horiz=F,beside=T,las=2,names.arg = c("Core","Contextual","Absent"),         
         ci.l =n.l,ci.u =n.u,
         col = brewer.pal(length(n),"RdBu")[c(1,3,2)],
         ylab="N° reactions",plot.ci=T,
         legend.text=F,ylim=c(0,4000))

#7: Terminate
write.table(rxns.core,"rxns_core.txt",sep = "\t",quote = F,row.names = F)
write.table(rxns.cntx,"rxns_cntx.txt",sep = "\t",quote = F,row.names =F)
write.table(rxns.absn,"rxns_absn.txt",sep = "\t",quote = F,row.names = F)