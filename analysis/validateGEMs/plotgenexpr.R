plotgenexpr <- function(genes,label){
  require(ggplot2)
  require(reshape2)
  E.top <- v$E[c(genes),]
  rownames(E.top) <- annotEntrez[rownames(E.top),"hgnc_symbol"]
  E.sel <- data.frame(t(E.top),
                      ctcl   = CTcl[!filter],
                      status = status[!filter],
                      ctclst = paste0(CTcl[!filter],"-",status[!filter]))
  E.sel$ctclst <- factor(E.sel$ctclst,
                         levels=c("ki-pa-re-co-normal","ki-pa-re-co-tumor",
                                  "lu-bl-br-hn-normal","lu-bl-br-hn-tumor",
                                  "ucec - ov-normal","ucec - ov-tumor"))
  E.m   <- melt(E.sel,measure.vars = 1:length(genes))  
  ggplot(E.m,aes(x=ctclst,y=value,fill=factor(ctclst))) +
    geom_boxplot() + xlab("") + ylab("Size-corrected log-cpm") +
    guides(fill=F)+ theme(axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))+
    facet_grid(~variable) + scale_fill_brewer(palette="Paired") +
    ggsave(paste0(label,"_cntxgenes_bwCL&vsN.pdf"),width=8,height=8)
}