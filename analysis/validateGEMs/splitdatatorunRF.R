splitdatatorunRF <- function(data,label,split=0.50,mtry=10,ntree=1000){
  library(randomForest)
  set.seed(1111)
  ind.test  <- c()
  ind.train <- c()
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
  return(list(rf=rf,data.test=data.test,label.test=label.test))
}