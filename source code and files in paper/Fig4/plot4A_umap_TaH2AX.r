library(reshape2)
library(Seurat)
library(ggplot2)
CS_scRNA <- readRDS("CS_mergereps.rds")
df <- as.data.frame(CS_scRNA@reductions$umap@cell.embeddings)
df$cluster <- as.factor(CS_scRNA@meta.data$seurat_clusters)
gene_df <- t(as.data.frame(GetAssayData(object = CS_scRNA,layer = "data")[c("TraesCS5A02G098300","TraesCS5B02G103600","TraesCS5D02G110600"),]))
df_out <- merge(gene_df,df,by=0,all=T)
df_melt1 <- reshape2::melt(df_out[,c("UMAP_1","UMAP_2","cluster","TraesCS5A02G098300","TraesCS5B02G103600","TraesCS5D02G110600")],id.vars=c("UMAP_1","UMAP_2","cluster"),measure.vars=c("TraesCS5A02G098300","TraesCS5B02G103600","TraesCS5D02G110600"))
df_melt1$valuenor <- ifelse(df_melt1$value>=0.5,1,df_melt1$value)
p <- ggplot(df_melt1,aes(UMAP_1,UMAP_2,color=valuenor))+geom_point(size=0.1)+
  scale_colour_gradientn(colours = c("#B7D2EA","white","#DD6434"))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "",x="UMAP 1",y="UMAP 2")
ggsave("plot4A_umap_of_singlecell_TaH2AX.pdf",p,width = 4.2,height = 3.5,dpi = 300)