performRFongenexpronlyT <- function(comparison,mtry=10,ntree=1000,randomize=F){
  #1: Choose genes to use to classify
  genes.sel <- cntx.genes.forcls
  
  #2: Get expression data and limit to selected genes
  load("~/Documents/Academia/Projects/MutLands/ProjectCodes/MutLands/validateResults/preprocessed_val.RData")
  E.T  <- y$E
  ct.T <- ctfactor
  genes.sel <- genes.sel[genes.sel %in% rownames(E.T)]
  data      <- t(E.T[genes.sel,])
  
  #3: Annotate sample for cluster and status
  CTcl.T <- as.character(ct.T)
  CTcl.T[CTcl.T=="lgg" | CTcl.T == "gbm"] <- "lgg - gbm"
  CTcl.T[CTcl.T=="ucec" | CTcl.T == "ov"] <- "ucec - ov"
  CTcl.T[CTcl.T=="kirc" | CTcl.T == "paad" | 
           CTcl.T == "read" | CTcl.T == "coad"] <- "ki-pa-re-co"
  CTcl.T[CTcl.T=="lusc" | CTcl.T == "luad" | 
           CTcl.T == "blca" | CTcl.T == "brca" |
           CTcl.T == "hnsc"] <- "lu-bl-br-hn"
  CTcl   <- factor(CTcl.T)
  
  #3: Split data into train and test and run RF
  source("splitdatatorunRF.R")
  ###First multilabel classification into clusters
  res       <- splitdatatorunRF(data,CTcl,mtry=mtry, ntree=ntree)
  ###Evaluate performance of multilabel split
  library(pROC)
  rf.pr     <- as.numeric(predict(res[["rf"]],type="response",
                                  newdata=res[["data.test"]]))
  rf.perf   <- multiclass.roc(res[["label.test"]],rf.pr)
  #plot the confusion matrix
  library(ggplot2)
  library(reshape2)
  source("pltoconfmat.R")
  confusion.train <- melt(res$rf$confusion[,-5])
  plotconfmat(confusion.train,"training")
  confusion.test  <- melt(res$rf$test$confusion[,-5])
  plotconfmat(confusion.test,"test")
  
  #5: Compute random classifiers
  if (randomize){
    library(doMC)
    library(foreach)
    set.seed(1234)
    registerDoMC(10)
    nsamples  <- 1000
    #Gene randomization
    classif.r <- foreach(b = 1:nsamples,.combine = cbind) %dopar% {
      # start <- proc.time()
      #Select random list of genes and get expres
      genes.sel.r <- sample(rownames(E.T),length(genes.sel),replace = T)  
      data.r      <- t(E.T[genes.sel.r,])
      #Split data into train and test and run RF
      res.r       <- splitdatatorunRF(data.r,CTcl,mtry=mtry, ntree=ntree)
      #Compute random performance on multiclass
      rf.pr.r     <- as.numeric(predict(res.r[["rf"]],type="response",
                                        newdata=res.r[["data.test"]]))
      rf.perf.r   <- multiclass.roc(res.r[["label.test"]],rf.pr.r)
      classif.r   <- as.numeric(rf.perf.r$auc)
      #     end <- proc.time()
      #     end-start
      classif.r} #Collect performance in a vector: multiclass + auc x cl
    rownames(classif.r) <- "Multiclass AUC"
    save("classif.r",file="random66geneclass_onlyT.rda")
    auc.notr <- as.numeric(rf.perf$auc)
    p.perm   <- sum(classif.r >= auc.notr)/nsamples
  
    #Sample randomization
    classif.s <- foreach(b = 1:nsamples,.combine = cbind) %dopar% {
      # start <- proc.time()
      #Select random list of genes and get expres
      CTcl.r <- sample(CTcl,length(CTcl),replace = F)  
      #Split data into train and test and run RF
      res.r       <- splitdatatorunRF(data,CTcl.r,mtry=mtry, ntree=ntree)
      #Compute random performance on multiclass
      rf.pr.r     <- as.numeric(predict(res.r[["rf"]],type="response",
                                        newdata=res.r[["data.test"]]))
      rf.perf.r   <- multiclass.roc(res.r[["label.test"]],rf.pr.r)
      classif.r   <- as.numeric(rf.perf.r$auc)
      #     end <- proc.time()
      #     end-start
      classif.r} #Collect performance in a vector: multiclass + auc x cl
    rownames(classif.s) <- "Multiclass AUC"
    save("classif.s",file="66generandomclass_onlyT.rda")
    auc.notr <- as.numeric(rf.perf$auc)
    p.perm.s <- sum(classif.s >= auc.notr)/nsamples
    
  #6: Plot density of random AUC
    classif.r.d <- melt(classif.r)
    ggplot(classif.r.d,aes(x=value,fill=Var1)) + 
      geom_density(alpha=0.3,bw="SJ") + theme_bw() +
      geom_vline(xintercept = as.numeric(rf.perf$auc),lwd=1.2) +
      guides(fill=F)+ xlab("Multiclass AUC") + ylab("Density") +
      ggsave("multiclassAUCranddensity_onlyT.pdf",width=14, height=8)
  }
  return(list(rf.res=res,rf.perf=rf.perf))
}