#Load list of good and bad models and all potential models
setwd("~/Documents/Academia/Projects/RxnLands/Codes/analyzeGEMs/QC/")
require(RColorBrewer)

#Load
raw <- read.delim("INITmodelSummary_goodModels.txt", sep="\t")
PSID.reconstructed.good <- rownames(raw)
raw <- read.delim("INITmodelSummary_badModels.txt", sep="\t")
PSID.reconstructed.bad  <- rownames(raw)
PSID.reconstructed.all  <- c(PSID.reconstructed.bad,PSID.reconstructed.good)
load("../../data/preprocessed.RData")
PSID.all                <- rownames(y$samples)

#Derive discarded PSID and define their features
PSID.notconverged    <- setdiff(PSID.all,PSID.reconstructed.all)
PSID.discarded       <- union(PSID.notconverged,PSID.reconstructed.bad)
filter               <- PSID.all%in%PSID.discarded
design.PSID.discarded <- data.frame(fGlobalCT[filter],
                                   fMuts[filter,c("APC","CASP8",
                                                       "CTNNB1","KEAP1",
                                                       "KRAS","IDH1","NFE2L2",
                                                       "NSD1","PTEN",
                                                       "RB1","STK11","TP53")])
cT.PSID.discarded      <- summary(design.PSID.discarded[,"fGlobalCT.filter."])
cT.PSID                <- table(fGlobalCT)
cT.PSID.discarded.perc <- cT.PSID.discarded/cT.PSID*100
cT.PSID.discarded.perc <- cT.PSID.discarded.perc[order(cT.PSID.discarded.perc)]
muts.PSID.discarded <- apply(fMuts[filter,c("APC","CASP8","CTNNB1",
                                                 "KEAP1","KRAS","IDH1","NFE2L2","NSD1","PTEN",
                                                 "RB1","STK11","TP53")]=="1",2,sum)
muts.PSID <- apply(fMuts[,c("APC","CASP8","CTNNB1",
                                 "KEAP1","KRAS","IDH1","NFE2L2","NSD1","PTEN",
                                 "RB1","STK11","TP53")]=="1",2,sum)
muts.PSID.discarded.perc <- muts.PSID.discarded/muts.PSID*100

#Plot
par(mfrow=c(1,1))
barplot(muts.PSID.discarded.perc,ylim=c(0,100),
        col=brewer.pal(length(names(muts.PSID.discarded.perc)),"Set3"),
        las=2,ylab="% of discarded PSID bearing the mutation")
bp <- barplot(cT.PSID.discarded.perc,ylim=c(0,100),names.arg = toupper(names(cT.PSID.discarded.perc)),
        col=brewer.pal(length(names(cT.PSID.discarded.perc)),"Set3"),
        las=2,ylab="% of discarded GEMs within GEMs of a cancer type")
text(bp,cT.PSID.discarded.perc+2,
     cT.PSID.discarded[order(cT.PSID.discarded/cT.PSID)])
