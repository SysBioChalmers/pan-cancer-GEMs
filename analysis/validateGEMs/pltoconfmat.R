plotconfmat <- function(confusion,set){
  ggplot(confusion,aes(x=Var1,y=Var2)) + 
    geom_tile(aes(fill=value)) + theme_bw() + 
    scale_fill_gradient2() + xlab("") + ylab("") + 
    geom_text(aes(label=value)) + ggtitle(paste0(set," set")) + 
    ggsave(paste0("confmat_",set,"_CTcl.pdf"),
           width = 14.7,height = 8.32)
  }