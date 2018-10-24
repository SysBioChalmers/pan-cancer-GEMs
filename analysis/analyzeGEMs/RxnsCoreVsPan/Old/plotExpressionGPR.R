plotExpressionGPR <- function(rxnID){
  #Check for data
  stopifnot(exists("E"))
  if(!exists("geneInclMat")){
    geneIncMatFile <- c("../../geneInclMat.txt")
    geneInclMat  <- read.delim(geneIncMatFile,header=T,sep="\t",row.names=1,allowEscapes=T)
    geneInclMat  <- geneInclMat[,-grep("X.1",colnames(geneInclMat))]
  }
  if(!exists("HMR3675.gpr")){
    HMR3675.gpr    <- read.delim("../../../docs/cHMR3765_rxnGeneMatrix.txt")
    rownames(HMR3675.gpr)<-HMR3675.gpr[,1]
    HMR3675.gpr  <- HMR3675.gpr[,-c(1,ncol(HMR3675.gpr))]
  }
  if(!exists("HMR3676.rxnTable")){
    HMR3676.rxnTable <- read.delim("../../../docs/cHMR3765_rxnTable.txt")
  }
  
  #Manipulate input data
  E.ensg            <- E
  colnames(E.ensg)  <- gsub(pattern="-",replacement=".",x=colnames(E.ensg))
  rownames(E.ensg)  <- annotEntrez[rownames(E),"ensembl_gene_id"]
  genes.ensg.common <- intersect(rownames(E.ensg),rownames(geneInclMat))
  geneInclMat.c     <- geneInclMat[genes.ensg.common,]
  E.ensg.c          <- E.ensg[genes.ensg.common,]
  E.ensg.c.zerneg   <- E.ensg.c
  E.ensg.c.zerneg[E.ensg.c<0] <- 0
  rxn.ex             <- as.character(rxnID)
  HMR3675.gpr.num    <- as.matrix(HMR3675.gpr)  
  
  #Retrieve gpr and plot expression
  rxn.ex.genes       <- colnames(HMR3675.gpr.num)[which(HMR3675.gpr.num[rxn.ex,]!=0)]
  rxn.ex.genes.hgnc  <- annotEntrez[match(rxn.ex.genes,annotEntrez[,"ensembl_gene_id"]),"hgnc_symbol"]
  rxn.ex.gene.sort   <- sort.int(E.ensg.c[rxn.ex.genes[1],],index.return = T)
  cols <- rainbow(length(rxn.ex.genes))
  barplot(E.ensg.c.zerneg[rxn.ex.genes,rxn.ex.gene.sort$ix],beside=F,las=2,
          col=cols,legend.text = rxn.ex.genes.hgnc,horiz=T,
          main=HMR3676.rxnTable[HMR3676.rxnTable[,"RxnID"]==rxn.ex,"Eq"],
          xlab="Size-corrected log-cpm")
}