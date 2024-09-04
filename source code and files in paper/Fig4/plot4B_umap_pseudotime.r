library(reshape2)
library(Seurat)
library(ggplot2)
pt <- read.table("pgDQ_monocle_pseudotime.xls",sep = "\t",stringsAsFactors = F,header = T)
p <- ggplot(pt,aes(UMAP_1,UMAP_2,color=Pseudotime))+geom_point(size=0.1)+
  scale_colour_gradientn(colours = c("#297AC2","white","#DD6434"))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "",x="UMAP 1",y="UMAP 2")
ggsave("umap_of_singlecell_Pseudotime.pdf",p,width = 4.2,height = 3.5,dpi = 300)