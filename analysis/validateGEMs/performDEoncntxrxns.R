performDEoncntxrxns <- function(comparison){
  #1: Get expression data and limit to selected genes
  load("~/Documents/Academia/Projects/MutLands/ProjectCodes/MutLands/validateResults/preprocessed_val.RData")
  y.T  <- y
  ct.T <- ctfactor
  load("~/Documents/Academia/Projects/MutLands/ProjectCodes/MutLands/validateResults/preprocessed_val_N.RData")
  y.N  <- y
  ct.N <- ctfactor

  #2: Creat sample factors for cluster and status
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
  status <- factor(c(rep("tumor",ncol(y.T)),rep("normal",ncol(y.N))))
  
  #3: Load gene expression in DGElist 
  require(limma)
  require(edgeR)
  y <- DGEList(counts=cbind(y.T$counts,y.N$counts))
  y <- calcNormFactors(y)
  
  #4: Model matrix and fit voom
  design <- model.matrix(~ 0 + CTcl : status)
  filter <- apply(design[,grepl("CTcllgg - gbm",colnames(design))]==1,1,any)  #Too few samples
    design <- design[!filter,!grepl("CTcllgg - gbm",colnames(design))]
    colnames(design) <- c("KPRC","LBBH","UO",
                          "KPRC_Tumor","LBBH_Tumor","UO_Tumor")
    y      <- y[,!filter]
  v      <- voom(y, design, plot = TRUE)
  fit    <- lmFit(v, design)
  contrast.matrix <- makeContrasts(contrasts=c("KPRC_Tumor - LBBH_Tumor",
                                               "LBBH_Tumor - UO_Tumor",
                                               "UO_Tumor - KPRC_Tumor",
                                               "KPRC - LBBH",
                                               "LBBH - UO",
                                               "UO - KPRC",
                                               "KPRC_Tumor - KPRC",
                                               "LBBH_Tumor - LBBH",
                                               "UO_Tumor - UO"),
                                   levels=design)
  fit2   <- contrasts.fit(fit, contrast.matrix)
  eb     <- eBayes(fit2)
  
  #5: Test for statistical significance
  coef2test <- colnames(contrast.matrix)
  p         <- as.matrix(eb$coef/eb$coef)
  for (coef2 in coef2test){
    tt <- topTable(eb,coef=coef2,number = Inf)
    p[,coef2] <- tt$adj.P.Val[match(rownames(p),rownames(tt))]}
  ranks     <- cbind(apply(p[,c("KPRC_Tumor - LBBH_Tumor",
                                "LBBH_Tumor - UO_Tumor",
                                "UO_Tumor - KPRC_Tumor")],2,rank),
                     apply(1/p[,c("KPRC - LBBH",
                                  "LBBH - UO",
                                  "UO - KPRC")],2,rank),
                     apply(p[,c("KPRC_Tumor - KPRC",
                                "LBBH_Tumor - LBBH",
                                "UO_Tumor - UO")],2,rank)
                     )
  
  #6: Choose genes to check
  genes.sel <- cntx.genes.forcls
  genes.sel <- genes.sel[genes.sel %in% rownames(v$E)]
  
#   ##More rel than random? Not really
#   deg.r     <- apply(matrix(1,nrow=length(genes.sel),ncol=10000),2,
#                     function(x){a<- sample(1:nrow(tt),length(genes.sel))
#                     x <- tt$adj.P.Val[a]
#                     return(x)})
#   prop.deg.r   <- colSums(deg.r<0.001)/length(genes.sel)
#   prop.deg.r.m <- melt(prop.deg.r)
#     sum(prop.deg.r >= prop.deg.sel)/100
#   ggplot(prop.deg.r.m,aes(x=value,fill="red")) + geom_density(alpha=0.3) +
#     theme_bw() + geom_vline(xintercept = prop.deg.sel)
  
  ##AWFUL CODING ALERT!
  #How many are deg in at least 1 tumor cluster comparison?
  #And how many are deg in at least 1 comparison with a matched normal tissue?
  signcontr <- matrix(nrow=length(rownames(p)),ncol=3,
                      dimnames=list(rownames(p),
                                    c("signinTCl","signvsN","notsigninN")))
  for (gene in rownames(p)){
    #Check if gene is de in at least 1 tumor cluster comparison
    signcl <- p[gene,c("KPRC_Tumor - LBBH_Tumor",
                       "LBBH_Tumor - UO_Tumor",
                       "UO_Tumor - KPRC_Tumor")]<0.001
    if (any(signcl)){
      pickcl   <- names(signcl)[signcl]
      #Is the difference between normal tissues relevant?
      notsignN <- p[gene,c("KPRC - LBBH",
                           "LBBH - UO",
                           "UO - KPRC")]>0.01
      if(any(notsignN)){
        pickN    <- names(notsignN)[notsignN]
        for (inN in pickN){
          cl <- paste0(paste0(unlist(strsplit(inN," - ")),"_Tumor")[1],
                 " - ",
                 paste0(unlist(strsplit(inN," - ")),"_Tumor")[2])
          if (cl %in% pickcl){
            signcontr[gene,1] <- cl
            if(is.na(signcontr[gene,2])){
              #If so, pick that cl comparison first and check if for each tumor cluster
              #in the comparison the gene is de vs. at least 1 match normal comparison
              #If not pick another significant comparison.
              #If not, leave NA          signcontr[gene,1] <- cl
              vsN.1 <- grep(substr(cl,1,6),colnames(p))[3]
              vsN.2 <- grep(substr(cl,14,18),colnames(p))[3]
              if (p[gene,vsN.1] < 0.001) {
                signcontr[gene,2] <- colnames(p)[vsN.1]
                signcontr[gene,3] <- inN
              } else if (p[gene,vsN.2] < 0.001) {
                signcontr[gene,2] <- colnames(p)[vsN.2]
                signcontr[gene,3] <- inN
              }
            }
          }
        }
      }
      #If not, pick a cl comparison and repeat
      if (is.na(signcontr[gene,1])){
        for (cl in pickcl){
          signcontr[gene,1] <- cl 
          if(is.na(signcontr[gene,2])){
            #If so, pick that cl comparison first and check if for each tumor cluster
            #in the comparison the gene is de vs. at least 1 match normal comparison
            #If not pick another significant comparison.
            #If not, leave NA          signcontr[gene,1] <- cl
            vsN.1 <- grep(substr(cl,1,6),colnames(p))[3]
            vsN.2 <- grep(substr(cl,14,18),colnames(p))[3]
            if (p[gene,vsN.1] < 0.001) {
              signcontr[gene,2] <- colnames(p)[vsN.1]
              signcontr[gene,3] <- inN
            } else if (p[gene,vsN.2] < 0.001) {
              signcontr[gene,2] <- colnames(p)[vsN.2]
              signcontr[gene,3] <- inN
            }
          }
        }
      }
    }
  }
  ndeg          <- colSums(!is.na(signcontr))/nrow(signcontr)
  signcontr.sel <- signcontr[genes.sel,]
  ndeg.sel      <- colSums(!is.na(signcontr.sel))
  
  #Detect interesting selected genes
  gene.sel.int  <- rownames(signcontr.sel)[apply(!is.na(signcontr.sel),1,all)]
  ranks.sel.int <- ranks[gene.sel.int,]
    tt <- topTable(eb,number = Inf)
  expr.sel      <- rownames(tt)[tt$AveExpr>0]
  top9          <- c("")
  for (ind in 1:3){
    inds       <- matrix(1:6,ncol=2)[ind,]
    rankvs     <- cbind(ranks.sel.int[,inds[1]],
                        ranks.sel.int[,inds[2]])
    toprankvs  <- rownames(rankvs[order(rowMeans(rankvs)),])
    top9[c(inds,6+ind)] <- toprankvs[toprankvs %in% expr.sel][1:3]
  }
  p.sel.top9 <- p[top9,]
  
  #6: Plot genes
  source("plotgenexpr.R")
  discussed  <- c("140838","200010","29881","219",
                  "217","224","501","223")
  plotgenexpr(discussed,"discussed")
  plotgenexpr(discussed,"top9")
}
