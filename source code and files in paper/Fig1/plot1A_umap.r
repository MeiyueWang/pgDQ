library(Seurat)
library(ggplot2)
CS_scRNA <- readRDS("CS_mergereps.rds")
dfABD <- data.frame(CS_scRNA$seurat_clusters,stringsAsFactors = F)
dfABD$CS_scRNA.seurat_clusters <- as.character(dfABD$CS_scRNA.seurat_clusters)
df <- as.data.frame(CS_scRNA@reductions$umap@cell.embeddings)
df$clusterABD <- dfABD[rownames(df),"CS_scRNA.seurat_clusters"]
df$clusterABD <- factor(df$clusterABD,levels = as.character(unique(sort(as.numeric(df$clusterABD)))))
df.color <- data.frame(no=c(0:19),color=c("#B3CAE8","#2077B5","#FF7F0E","#FFBD7E","#37A537","#98DF8B","#D52728","#FF9E9C","#9467BC","#C5B1D5","#905B52","#C59B96","#E276C1","#F7B9D4","#7E7E7E","#C6C6C6","#C5C53D","#CACA52","#DBDB8D","#939B6F"))
p <- ggplot(df,aes(UMAP_1,UMAP_2,color=clusterABD))+geom_point(size=0.5)+
  scale_color_manual(values = df.color$color)+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "",x="UMAP 1",y="UMAP 2")
ggsave("singlecell_UMAP.pdf",p,height = 4,width = 5,units = "in",dpi = 300)