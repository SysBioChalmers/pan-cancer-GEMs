writeHistogramFile <- function(mutGeneMatrix.met.sign.sort,genes.maxSize,
                               filename,scale=1,dir=T){
  #If there is a hist file, replace it.
  if (file.exists(filename)){file.remove(filename)}
  file.hm <- file(filename,"at")
  #Go through each gene
  for (gene in rownames(mutGeneMatrix.met.sign.sort)){
    #Count n° of upreg muts
    upcount <- sum(mutGeneMatrix.met.sign.sort[gene,]>0)
    #Count n° of dnreg muts
    dncount <- sum(mutGeneMatrix.met.sign.sort[gene,]<0)
    #Scale width of hist bar to the n° of muts
    scaleFactor <- scale*sum(mutGeneMatrix.met.sign.sort[gene,]!=0)
    #Use directions?
    if (dir){
    #Are there any muts associated?
    if (sum(mutGeneMatrix.met.sign.sort[gene,]!=0)>1){
      tick        <- floor(genes.maxSize*(scaleFactor + abs(scale-1))/2)
      string <- paste(gene,0,tick,upcount,paste0("id=","up"),sep="\t")
      writeLines(string,file.hm)
      string <- paste(gene,tick+1,genes.maxSize*(scaleFactor + abs(scale-1)),dncount,paste0("id=","dn"),sep="\t")
      writeLines(string,file.hm)} else {
        if (dncount==1){
        string <- paste(gene,0,genes.maxSize*(scaleFactor + abs(scale-1)),dncount,paste0("id=","dn"),sep="\t")
        writeLines(string,file.hm)} else {
          string <- paste(gene,0,genes.maxSize*(scaleFactor + abs(scale-1)),upcount,paste0("id=","up"),sep="\t")
          writeLines(string,file.hm)        
        }
      }
  } else {
  #Don't use directions
    nregmuts <- sum(mutGeneMatrix.met.sign.sort[gene,]!=0)
    string <- paste(gene,0,genes.maxSize,nregmuts+1,
                    paste0("fill_color=blues-9-seq-",nregmuts+1),sep="\t")
    writeLines(string,file.hm)
  
  }}
  close(file.hm)}