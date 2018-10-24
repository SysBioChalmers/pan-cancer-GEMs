pcaElInclMatrix <-function(rxnInclMat,loads=F){
  ### Cluster gene in models
  require(ade4)
  require(RColorBrewer)
  E.pca  <- dudi.pca(t(rxnInclMat),center=T,scale=FALSE,scan=FALSE)
  totVar <- sum(E.pca$eig)
  var1   <- E.pca$eig[1]/totVar*100
  var2   <- E.pca$eig[2]/totVar*100
  if(loads){par(mfrow=c(2,1))}
  s.class(E.pca$li,fGlobalCT.m,cpoint=1,
          col=c(brewer.pal(nCancerTypes,"Set3"),"light blue"),
          label=toupper(levels(fGlobalCT.m)))
  text(10,0.5,paste0(round(var1,1),"%"))
  text(-0.5,5,paste0(round(var2,1),"%"),srt=90)
  if(loads){
    s.arrow(E.pca$c1,clabel=1)
    par(mfrow=c(1,1))}    
}