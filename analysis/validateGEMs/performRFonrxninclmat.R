performRFonrxninclmat <- function(comparison,mtry=10,ntree=1000){
  library(pROC)
  library(randomForest)
  #1: Choose rxns to use to classify
    if (comparison=="TvsN"){
      rxns.sel    <- core.rxns.onlyT
      data        <- rbind(t(rxnInclMat[HMR2.rxnTable[rxns.sel,"Eq"],]),
                         t(rxnInclMat_N[rxns.sel,]))}
    if (comparison=="Tclus"){
      rxns.sel    <- cntx.rxns.forcls
      data        <- t(rxnInclMat[HMR2.rxnTable[rxns.sel,"Eq"],])}
  
    #2: Assign label
    if (comparison=="Tclus"){
      label <- as.character(fGlobalCT.m)
      label[label=="lgg" | label == "gbm"] <- "lgg - gbm"
      label[label=="ucec" | label == "ov"] <- "ucec - ov"
      label[label=="kirc" | label == "paad" | 
              label == "read" | label == "coad"] <- "ki-pa-re-co"
      label[label=="lusc" | label == "luad" | 
              label == "blca" | label == "brca" |
              label == "hnsc"] <- "lu-bl-br-hn"}
    if (comparison=="TvsN"){
      label <- c(rep("tumor",ncol(rxnInclMat)),
                 rep("normal",ncol(rxnInclMat_N)))}
    label <- as.factor(label)
  
    #3: Split data into train and test: find indexes
    set.seed(1111)
    ind.test  <- c()
    ind.train <- c()
    split     <- 0.50
    for (clust in levels(label)){
      ind.clinlab <- which(label %in% clust)
      ind.test.cl    <- sample(ind.clinlab,length(ind.clinlab)*split,replace = F)
      ind.train.cl   <- setdiff(ind.clinlab,ind.test.cl)
      ind.test  <- c(ind.test,ind.test.cl)
      ind.train <- c(ind.train,ind.train.cl)
    }
    data.train  <- data[ind.train,]
    data.test   <- data[ind.test,]
    label.train <- label[ind.train]
    label.test  <- label[ind.test]
    rf          <- randomForest(data.train,label.train,
                                data.test,label.test,
                                mtry=mtry, ntree=ntree,
                                keep.forest=TRUE, importance=TRUE)
    rf.pr       <- as.numeric(predict(rf,type="response",
                                      newdata=data.test))
    multi.roc   <- multiclass.roc(label.test,rf.pr)
    
    res     <- list(rf=rf,multi.roc=multi.roc)
    return(res)
}