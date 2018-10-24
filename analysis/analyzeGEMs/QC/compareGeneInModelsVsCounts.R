### Compare gene read counts with gene inclusion in PSID models
rm(list=ls())
require(ggplot2)
require(reshape2)
setwd("~/Documents/Academia/Projects/RxnLands/Codes/analyzeGEMs/QC/")

# Load data
load("../../data/annotation.RData")
load("../../data/preprocessed.RData")
geneInclMat <- read.delim("../geneInclMat.txt",header=T,sep="\t",row.names=1,allowEscapes=T)
geneInclMat <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
model.summary <- read.delim("../INITmodelSummary.txt",row.names="PSID")

# Normalize gene inclusion matrix with read count matrix
lib.size          <- with(y$samples, lib.size * norm.factors)
E                 <- t(log2(t(y$counts + 0.5)/(lib.size + 1) * 1e+06))
E.ensg            <- E
colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat))
geneInclMat.c     <- geneInclMat[genes.ensg.common,]
E.ensg.c          <- E.ensg[genes.ensg.common,colnames(geneInclMat.c)]

# Join info about expr value and inclusion per sample
m <- melt(t(E.ensg.c))
n <- melt(t(geneInclMat.c))
  ##Check if identical variable
  identical(n$Var1,m$Var1)
  identical(n$Var2,m$Var2)
o <- data.frame("PSID"=m$Var1,"Gene"=m$Var2,"E"=m$value,"Included"=n$value)

#Density plots
##Generic
ggplot(o,aes(x=E,fill=factor(Included))) + geom_density() +
    geom_density(alpha=.5,size=0) +  ylab("Density") +
    scale_fill_manual(values=c(brewer.pal(9,"Reds")[5],
                               brewer.pal(9,"Blues")[5])) +
    xlab("Size-adjusted log-cpm") + theme_bw() + 
    theme(legend.position="none") + ylim(0,0.3) + xlim(-5,15) +
    ggsave(filename = "Eleveldens_Excl_allmodels.pdf",
         width=10,height=6)

##Per sample
set.seed(1111)
PSID.random <- sample(unique(o$PSID),9,F)
ggplot() + 
  geom_density(data=subset(o,Included==1 & PSID %in% PSID.random),
               aes(x=E,fill=factor(PSID)),alpha=.3,size=0) +  
  scale_fill_brewer(palette="Blues") +
  ylab("Density") + xlab("Size-adjusted log-cpm") + theme_bw() + 
  theme(legend.position="none") + ylim(0,0.3) + xlim(-5,15) +
  ggsave(filename = "Eleveldens_Incl_xmodels.pdf",
         width=10,height=6)

ggplot() + 
  geom_density(data=subset(o,Included==0 & PSID %in% PSID.random),
               aes(x=E,fill=factor(PSID)),alpha=.3,size=0) +  
  scale_fill_brewer(palette="Reds") +
  ylab("Density") + xlab("Size-adjusted log-cpm") + theme_bw() + 
  theme(legend.position="none") + ylim(0,0.3) + xlim(-5,15) +
  ggsave(filename = "Eleveldens_Excl_xmodels.pdf",
         width=10,height=6)

# Boxplots for gene In vs counts for 6 random models
random          <- sample(x=1:dim(geneInclMat)[2],6,F)
PSID.example    <- colnames(geneInclMat)[random]
par(mfrow=c(2,3))
for (PSID in PSID.example){
  boxplot(E.ensg.c[,PSID]~geneInclMat.c[,PSID],col=c("Dark Red","Dark Blue"),
    names=c("ExcludedGenes","IncludedGenes"),ylab="Normalized log-cpm",
    main=PSID)}

# Empirical cumulative distributions for 4 random models
random          <- sample(x=1:dim(geneInclMat)[2],4,F)
PSID.example    <- colnames(geneInclMat)[random]
par(mfrow=c(4,1))
for (PSID in PSID.example){
  included <-  as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==1]
  excluded <-  as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==0]
  included_ecdf <- ecdf(included)
  excluded_ecdf <- ecdf(excluded)
  plot(included_ecdf,col="Dark Blue",xlab="Normalized log-cpm",
       ylab="Cumulative fraction",main=PSID)
  lines(excluded_ecdf,col="Dark Red")
  abline(h=0.5,col="Red")}

# Boxplots and ecdf for gene In vs counts in all models
included <- c()
excluded <- c()
k <- 0
for (PSID in colnames(geneInclMat)){
  k <- k+1
  print(k)
  included <- c(included,
                as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==1])
  excluded <- c(excluded,
                as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==0])
}
status  <- as.factor(c(rep(0,length(excluded)),rep(1,length(included))))
logcpms <- c(excluded,included)

### Find value that best separates included from excluded genes based on log-cpms
fit  <- rpart(status ~ logcpms)
pfit <- prune(fit, cp=fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
sep  <- 1.609
par(mfrow=c(1,1))
plot(pfit, uniform=TRUE, main="Pruned Classification Tree for Gene Inclusion")
  text(pfit, use.n=TRUE, all=TRUE, cex=.8)

boxplot(logcpms ~ status,
        col=c("Dark Red","Dark Blue"),
        names=c("ExcludedGenes","IncludedGenes"),ylab="Normalized log-cpm",
        main="All Models")
  abline(h=sep,col="black",lwd=2)

included_ecdf <- ecdf(included)
excluded_ecdf <- ecdf(excluded)
plot(included_ecdf,col="Dark Blue",xlab="Normalized log-cpm",
     ylab="Cumulative fraction",main="All Models")
lines(excluded_ecdf,col="Dark Red")
abline(h=0.5,col="Red")

# Correlate distance of cumulative distributions between IN and OUT in each model w accuracy stats
D <- matrix(0,nrow=dim(geneInclMat)[2],ncol=1,dimnames=list(colnames(geneInclMat),"KSstatistic"))
for (PSID in colnames(geneInclMat)){
  included <-  as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==1]
  excluded <-  as.vector(as.matrix(E.ensg.c[,PSID]))[geneInclMat.c[,PSID]==0]
  kstest <- ks.test(excluded,included,alternative="greater") #It saturates before cos it containes more neg values than "included"
  D[PSID,] <- kstest$statistic}
model.stats <- model.summary[,c("MCC","RxnAddedForTask...")]
rownames(model.stats) <- gsub(pattern="-",replacement=".",x=rownames(model.stats))
model.stats.f <- model.stats[match(rownames(D),rownames(model.stats)),]
par(mfrow=c(1,2))
  plot(D[,"KSstatistic"],model.stats.f[,"MCC"],
       xlab="Distance in gene expression between included vs excluded genes in a model",
       ylab="Accuracy of model reconstruction")
  text(x=(max(D[,"KSstatistic"])-0.1*max(D[,"KSstatistic"])),
       y=(min(model.stats.f[,"MCC"])+0.01*min(model.stats.f[,"MCC"])),
       labels=cor(D[,"KSstatistic"],model.stats.f[,"MCC"]))
  plot(D[,"KSstatistic"],model.stats.f[,"RxnAddedForTask..."],
     xlab="Distance in gene expression between included vs excluded genes in a model",
     ylab="N° rxns added for tasks")
  text(x=(max(D[,"KSstatistic"])-0.1*max(D[,"KSstatistic"])),
     y=min(model.stats.f[,"RxnAddedForTask..."])+0.1*min(model.stats.f[,"RxnAddedForTask..."]),
       labels=cor(D[,"KSstatistic"],model.stats.f[,"RxnAddedForTask..."]))  