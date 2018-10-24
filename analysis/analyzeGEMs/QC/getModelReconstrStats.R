###Filter reconstructed models upon quality
case <- "tumor"
#Get stats
if (case =="tumor"){
  model.summary <- read.delim("../INITmodelSummary.txt",row.names="PSID")
  suffix <- ""
} else {
  model.summary <- read.delim("../INITmodelSummary_N.txt",row.names="PSID")
  suffix <- "_N"
}
  par(mfrow=c(1,2))
  boxplot(model.summary[,c("ACC","MCC")],col=c("Blue","Dark Blue"),ylim=c(0,1))
  boxplot(model.summary[,"RxnAddedForTask..."],ylab="% of reaction in model added for tasks",col="grey",ylim=c(0,5))

#Draw density distributions
require(ggplot2)
require(reshape2)
data <- melt(model.summary,measure.vars = c("ACC","MCC","RxnAddedForTask..."))
ggplot(data,aes(x=value)) + geom_density() + facet_grid(~variable) +
  theme_bw() + xlim(0,1) + ylab("Density")
                  
#Fish good models
criterion1 <- model.summary[,"ACC"]>0.90
criterion2 <- model.summary[,"MCC"]>0.80
criterion3 <- model.summary[,"RxnAddedForTask..."] < 1
models.good.number <- sum(criterion1 & criterion2 & criterion3)
write.table(model.summary[criterion1 & criterion2 & criterion3,],
            file=paste0("INITmodelSummary_goodModels",suffix),quote=F,sep="\t")
write.table(model.summary[!(criterion1 & criterion2 & criterion3),],
            file=paste0("INITmodelSummary_badModels",suffix),quote=F,sep="\t")

#Print good models
models.good.names.MAT <- model.summary[(criterion1 & criterion2 & criterion3),
                                      c("ModelNameMAT")]
models.good.names.XML <- model.summary[(criterion1 & criterion2 & criterion3),
                                      c("ModelNameXML")]
write.table(models.good.names.MAT,file=paste0("listGoodReconstructedModelMATnames",suffix),
            quote=F,sep="\t",row.names = F,col.names = F)
write.table(models.good.names.XML,file=paste0("listGoodReconstructedModelXMLnames",suffix),
            quote=F,sep="\t",row.names = F,col.names = F)

#Print bad models
models.bad.names.MAT <- model.summary[!(criterion1 & criterion2 & criterion3),
                                    c("ModelNameMAT")]
models.bad.names.XML <- model.summary[!(criterion1 & criterion2 & criterion3),
                                    c("ModelNameXML")]
write.table(models.bad.names.MAT,file=paste0("listBadReconstructedModelMATnames",suffix),
            quote=F,sep="\t",row.names = F,col.names = F)
write.table(models.bad.names.XML,file=paste0("listBadReconstructedModelXMLnames",suffix),
            quote=F,sep="\t",row.names = F,col.names = F)