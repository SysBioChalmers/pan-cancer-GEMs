performRFongenexpr <- function(comparison,mtry=10,ntree=1000){
  #1: Choose genes to use to classify
  if (comparison=="TvsN"){genes.sel  <- core.genes.onlyT}
  if (comparison=="Tclus"){genes.sel <- cntx.genes.forcls}
  
  #2: Get expression data and limit to selected genes
  load("~/Documents/Academia/Projects/MutLands/ProjectCodes/MutLands/validateResults/preprocessed_val.RData")
  E.T  <- y$E
  ct.T <- ctfactor
  load("~/Documents/Academia/Projects/MutLands/ProjectCodes/MutLands/validateResults/preprocessed_val_N.RData")
  E.N  <- y$E
  ct.N <- ctfactor
  genes.sel <- genes.sel[genes.sel %in% rownames(E.T)]
  data      <- rbind(t(E.T[genes.sel,]),t(E.N[genes.sel,]))
  
  #3: Annotate sample for cluster and status
  CTcl.T <- as.character(ct.T)
    CTcl.T[CTcl.T=="lgg" | CTcl.T == "gbm"] <- "lgg - gbm"
    CTcl.T[CTcl.T=="ucec" | CTcl.T == "ov"] <- "ucec - ov"
    CTcl.T[CTcl.T=="kirc" | CTcl.T == "paad" | 
            CTcl.T == "read" | CTcl.T == "coad"] <- "ki-pa-re-co"
    CTcl.T[CTcl.T=="lusc" | CTcl.T == "luad" | 
            CTcl.T == "blca" | CTcl.T == "brca" |
            CTcl.T == "hnsc"] <- "lu-bl-br-hn"
  CTcl.N <- as.character(ct.N)
    CTcl.N[CTcl.N=="lgg" | CTcl.N == "gbm"] <- "lgg - gbm"
    CTcl.N[CTcl.N=="ucec" | CTcl.N == "ov"] <- "ucec - ov"
    CTcl.N[CTcl.N=="kirc" | CTcl.N == "paad" | 
              CTcl.N == "read" | CTcl.N == "coad"] <- "ki-pa-re-co"
    CTcl.N[CTcl.N=="lusc" | CTcl.N == "luad" | 
              CTcl.N == "blca" | CTcl.N == "brca" |
              CTcl.N == "hnsc"] <- "lu-bl-br-hn"    
  CTcl   <- factor(c(CTcl.T,CTcl.N))
  status <- factor(c(rep("tumor",ncol(E.T)),rep("normal",ncol(E.N))))
  
  #3: Split data into train and test and run RF
  source("splitdatatorunRF.R")
  ###First multilabel classification into clusters
  res  <- splitdatatorunRF(data,CTcl,mtry=mtry, ntree=ntree)
  ###Evaluate performance of multilabel split
  library(pROC)
  rf.pr     <- as.numeric(predict(res[["rf"]],type="response",
                                  newdata=res[["data.test"]]))
  rf.perf   <- multiclass.roc(res[["label.test"]],rf.pr)
  #plot the confusion matrix
  library(ggplot2)
  library(reshape2)
  confusion.train <- melt(res$rf$confusion[,-5])
  plotconfmat(confusion.train,"training")
  confusion.test  <- melt(res$rf$test$confusion[,-5])
  plotconfmat(confusion.test,"test")
  ###Then tumor-normal split within cluster
  pdf("ROCforTvsN_withinCL.pdf",width = 14.7,height = 4)
  par(mfrow=c(1,4))
  for (cluster in levels(CTcl)){
    res.cl  <- splitdatatorunRF(data[CTcl==cluster,],
                                status[CTcl==cluster],
                                mtry=mtry, ntree=ntree)
  ###Evaluate performance of binary split
    library(ROCR)
    rf.pr   <- predict(res.cl[["rf"]],type="prob",
                       newdata=res.cl[["data.test"]])[,2]
    rf.pred <- prediction(rf.pr, res.cl[["label.test"]])
    #performance in terms of true and false positive rates
    rf.perf = performance(rf.pred,"tpr","fpr")
    #compute area under curve
    auc <- performance(rf.pred,"auc")
    auc <- unlist(slot(auc, "y.values"))
    #plot the curve
    plot(rf.perf,col=2,lwd=2,main = cluster)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.1,paste0("AUC = ",round(auc,3)))
  }
  dev.off()
  
  #5: Compute random classifiers
  library(doMC)
  library(foreach)
  set.seed(1234)
  registerDoMC(10)
  nsamples  <- 1000
  rf.perf.r <- foreach(b = 1:nsamples,.combine = cbind) %dopar% {
    # start <- proc.time()
    #Select random list of genes and get expres
    genes.sel.r <- sample(rownames(E.T),length(genes.sel),replace = T)  
    data.r      <- rbind(t(E.T[genes.sel.r,]),t(E.N[genes.sel.r,]))
    #Split data into train and test and run RF
    res.r       <- splitdatatorunRF(data.r,CTcl,mtry=mtry, ntree=ntree)
    #Compute random performance on multiclass
    rf.pr.r     <- as.numeric(predict(res.r[["rf"]],type="response",
                                      newdata=res.r[["data.test"]]))
    rf.perf.r   <- multiclass.roc(res.r[["label.test"]],rf.pr.r)
    classif.r   <- as.numeric(rf.perf.r$auc)
    #Compure random performance on single CT
    for (cluster in levels(CTcl)){
      res.cl.r  <- splitdatatorunRF(data.r[CTcl==cluster,],
                                    status[CTcl==cluster],
                                    mtry=mtry, ntree=ntree)
      ###Evaluate performance of rndom binary split
      rf.pr.cl.r   <- predict(res.cl.r[["rf"]],type="prob",
                         newdata=res.cl.r[["data.test"]])[,2]
      rf.pred.cl.r <- prediction(rf.pr.cl.r, res.cl.r[["label.test"]])
      #compute area under curve of random classifier
      auc.cl.r <- performance(rf.pred.cl.r,"auc")
      auc.cl.r <- unlist(slot(auc.cl.r, "y.values"))
      classif.r<- c(classif.r,auc.cl.r)}
#     end <- proc.time()
#     end-start
    classif.r} #Collect performance in a vector: multiclass + auc x cl
  rownames(rf.perf.r) <- c("Multiclass AUC",
                           paste("AUC",levels(CTcl)))
  
  if (comparison=="TvsN"){
    pdf(file=paste0("ROCfor",comparison,"+random.pdf"),width = 6,height = 6)
    for (rf.perf.rx in rf.perf.r){
      plot(rf.perf.rx,col="gray")
      par(new=T)}
    plot(rf.perf,col=2,lwd=2)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.1,paste0("AUC = ",round(auc,3)))
    dev.off()}
}