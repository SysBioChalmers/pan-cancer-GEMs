writeHeatmapFile <- function(mutGeneMatrix.met.sign.sort,genes.maxSize,KEGG,filename,scale=1){
  if (file.exists(filename)){file.remove(filename)}
  file.hm <- file(filename,"at")
  for (gene in rownames(mutGeneMatrix.met.sign.sort)){
    if (gene %in% KEGG){heat=1} else {heat=0}
    scaleFactor <- scale*sum(mutGeneMatrix.met.sign.sort[gene,]!=0)
    string <- paste(gene,0,genes.maxSize*(scaleFactor + abs(scale-1)),heat,paste0("id=",gene),sep="\t")
    writeLines(string,file.hm)}
  close(file.hm)}