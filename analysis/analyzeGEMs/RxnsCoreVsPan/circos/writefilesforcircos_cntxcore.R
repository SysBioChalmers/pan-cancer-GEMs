###Color iPath map
rm(list=ls())
setwd("~/Box Sync/3rd Semester (May 2013 - November 2013)/Projects/RxnLands/Codes/analyzeGEMs/RxnsCoreVsPan/circos/")
require(RColorBrewer)

#Load
cntx.rxns   <- read.delim("../rxns_cntx.txt",stringsAsFactors=F)
coreiso.rxns   <- read.delim("../rxns_core_iso.txt",stringsAsFactors=F,header = F)
corepan.rxns   <- read.delim("../rxns_core_pan.txt",stringsAsFactors=F,header = F)
coreconn.rxns  <- read.delim("../rxns_core_conn.txt",stringsAsFactors=F,header = F)

#Write first karyo - core
filename.k <- paste0("data/karyotype_core.txt")
if (file.exists(filename.k)){file.remove(filename.k)}
file.karyo <- file(filename.k,"a+")
string <- paste("chr","-","corepan","pan",0,nrow(corepan.rxns),"ylorbr-3-seq-3",
                  sep="\t")
writeLines(string,file.karyo)
string <- paste("chr","-","coreiso","iso",0,nrow(coreiso.rxns),"ylorbr-3-seq-2",
                sep="\t")
writeLines(string,file.karyo)
string <- paste("chr","-","coreconn","conn",0,nrow(coreconn.rxns),"ylorbr-3-seq-1",
                sep="\t")
writeLines(string,file.karyo)
string <- paste("chr","-","cntx","cntx",0,nrow(cntx.rxns),"rdbu-3-div-3",
                sep="\t")
writeLines(string,file.karyo)
close(file.karyo)

