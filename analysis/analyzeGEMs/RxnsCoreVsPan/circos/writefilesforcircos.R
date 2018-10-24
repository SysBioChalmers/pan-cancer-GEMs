#Write circos file based on mut-met gene association
rm(list=ls())
require(limma)
require(gplots)
require(piano)
require(RColorBrewer)
require(synapseClient)
synapseLogin("gatto@chalmers.se","forre$ter1867")

#Load
annotation.file <- synGet("syn3168036")
load(annotation.file@filePath)
mutGeneMat.met.file    <- synGet("syn3192396")
mutGeneMatrix.met.sign <- read.delim(mutGeneMat.met.file@filePath)
metgenes.file   <- synGet("syn3192387")
metgenes        <- readLines(metgenes.file@filePath)
araxannot.file  <- synGet("syn3240888")
arax.annot      <- read.delim(araxannot.file@filePath,
                              sep="\t",stringsAsFactors=F)
writeheatm  <- synGet("syn3308147")
source(writeheatm@filePath)
writehist  <- synGet("syn3308156")
source(writehist@filePath)

#Parameters
convergence.cutoff <- 1 #Maximum n° of converging mutation to the gene to be shown
simplified         <- F
suffix             <- ifelse(simplified,"_simple","")

#1: Sort associations by number of converging mutations
mut.sort   <- sort.int(apply(mutGeneMatrix.met.sign!=0,2,sum),decreasing=T,index.return=T)
gene.sort  <- sort.int(apply(mutGeneMatrix.met.sign!=0,1,sum),decreasing=T,index.return=T)
mutGeneMatrix.met.sign.sort <- mutGeneMatrix.met.sign[gene.sort$ix,mut.sort$ix]
mutFactors <- paste0("metadata",colnames(mutGeneMatrix.met.sign.sort)) #This is done for retro compatibility
colnames(mutGeneMatrix.met.sign.sort) <- mutFactors

#2: Filter convergence lower than cutoff
mutGeneMatrix.met.sign.sort <- mutGeneMatrix.met.sign.sort[
  apply(mutGeneMatrix.met.sign.sort!=0,1,sum)>=convergence.cutoff,]

#3: Choose drawing settings
genes.maxSize <- 1000
muts.maxSize  <- apply(mutGeneMatrix.met.sign.sort!=0,2,sum)
ribbonScale   <- 0.5
col           <- paste0("set3-12-qual-",1:length(mutFactors))
names(col)    <- paste0("metadata",colnames(mutGeneMatrix.met.sign))#[sample(1:12,12,F)] #Assign colors randomly

#4: Write karyotype file
filename.k <- paste0("../circos/data/karyotype_mutations_allMets",suffix,".txt")
if (file.exists(filename.k)){file.remove(filename.k)}
file.karyo <- file(filename.k,"a+")
for (mutFactor in mutFactors){
  string <- paste("chr","-",mutFactor,substring(mutFactor,9),0,muts.maxSize[mutFactor],col[mutFactor],sep="\t")
  writeLines(string,file.karyo)}
for (gene in rownames(mutGeneMatrix.met.sign.sort)){
  if (gene!="3417"){#Skip IDH1
    gene.hgnc <- annotEntrez[gene,"hgnc_symbol"] #Get gene hgnc_symbol
    if (gene.hgnc==""){gene.hgnc<-annotEntrez[gene,"ensembl_gene_id"]} #If unavailable, use ensembl
    if (simplified){ #Same size for each gene, change name to arax genes
      string <- paste("chr","-",gene,gene.hgnc,0,genes.maxSize,"grey",sep="\t")
    } else { #Size proportional to convergence
      string <- paste("chr","-",gene,gene.hgnc,0,genes.maxSize*sum(mutGeneMatrix.met.sign.sort[gene,]!=0),"grey",sep="\t") #Scale according to convergence
    }
    writeLines(string,file.karyo)}}
close(file.karyo)

#5: Write link file
filename.l <- paste0("../circos/data/links_mutations_allMets",suffix,".txt")
if (file.exists(filename.l)){file.remove(filename.l)}
file.links <- file(filename.l,"at")
gene.count     <- apply(mutGeneMatrix.met.sign.sort!=0,1,sum)
gene.count.max <- gene.count
for (mutFactor in mutFactors){
  gene.muts <- rownames(mutGeneMatrix.met.sign.sort[mutGeneMatrix.met.sign.sort[,mutFactor]!=0,])
  k <- 0
  for (gene.mut in gene.muts){
    if (gene.mut!="3417"){
      fc <- 1 #Same thickness for every link
      if (simplified){#Cram all links to the gene ideogram
        string <- paste(mutFactor,k,k+1,gene.mut,
                        round(genes.maxSize*(gene.count[gene.mut]-1)/gene.count.max[gene.mut]),
                        round(genes.maxSize*(gene.count[gene.mut])/gene.count.max[gene.mut]),
                        paste0("color=",col[mutFactor],",thickness=",1+abs(fc),"p"),
                        sep="\t")        
      } else {#Span links along the gene ideogram
        string <- paste(mutFactor,k,k+1,gene.mut,
                        round(genes.maxSize*(gene.count[gene.mut]-1)),#/gene.count.max[gene.mut]),
                        round(genes.maxSize*(gene.count[gene.mut])),#/gene.count.max[gene.mut]),
                        paste0("color=",col[mutFactor],",thickness=",1+abs(fc),"p"),
                        sep="\t")}
      writeLines(string,file.links)
      k <- k+1 
      gene.count[gene.mut]<-gene.count[gene.mut]-1} else {#This is only for IDH1
        fc <- mutGeneMatrix.met.sign.sort[gene.mut,mutFactor]
        string <- paste(mutFactor,k,k+1,"metadataIDH1",
                        round(genes.maxSize*(gene.count[gene.mut]-1)),#/gene.count.max[gene.mut]),
                        round(genes.maxSize*(gene.count[gene.mut])),#/gene.count.max[gene.mut]),
                        paste0("color=",col[mutFactor],",thickness=",1+abs(fc),"p"),
                        sep="\t")
        writeLines(string,file.links)
        k <- k+1
        gene.count[gene.mut]<-gene.count[gene.mut]-1}}}
close(file.links)

#6: Write heatmap file
arax.genes <- rownames(arax.annot)
writeHeatmapFile(mutGeneMatrix.met.sign.sort,genes.maxSize,scale=ifelse(simplified,1,genes.maxSize),
                 arax.genes,paste0("../circos/data/heatmap_mutations_AraX",suffix,".txt"))

#7: Write histogram file
writeHistogramFile(mutGeneMatrix.met.sign.sort,genes.maxSize,scale=1,dir=T,
                   filename = "../circos/data/histogram_mutations_regdirec.txt")
writeHistogramFile(mutGeneMatrix.met.sign.sort,genes.maxSize,scale=1,dir=F,
                   filename = "../circos/data/histogram_mutations_regcnt.txt")

#8: Terminate
circos.data <- synGet("syn3192419")
karyo <- File(path = "../circos/data/karyotype_mutations_allMets.txt", parentId =circos.data$properties$id)
synSetAnnotations(karyo)<-list(description="A karyotype
                                          text file compatible with
                                          circos software for mutation
                                          and associated regulated
                                          metabolic genes")
karyo <- synStore(karyo)
links <- File(path = "../circos/data/links_mutations_allMets.txt", parentId =circos.data$properties$id)
synSetAnnotations(links)<-list(description="A link
                                          text file compatible with
                                          circos software for mutation
                                          and associated regulated
                                          metabolic genes.")
links <- synStore(links)
heatm <- File(path = "../circos/data/heatmap_mutations_AraX.txt", parentId =circos.data$properties$id)
synSetAnnotations(heatm)<-list(description="A heatmap
                               text file compatible with
                               circos software for metabolic genes
                               belonging to AraX")
heatm <- synStore(heatm)
hist <- File(path = "../circos/data/histogram_mutations_regdirec.txt", parentId =circos.data$properties$id)
synSetAnnotations(hist)<-list(description="A histogram
                                          text file compatible with
                                          circos software for mutation
                                          and associated regulated
                                          metabolic genes indicating
                                          n° up/dn.regulating mutations.")
hist <- synStore(hist)