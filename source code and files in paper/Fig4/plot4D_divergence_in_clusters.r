library(ggplot2)
pt <- read.table("pseudotime_divergence_in_cells.xls",sep = "\t",stringsAsFactors = F,header = T)
df.color <- data.frame(no=c(0:19),color=c("#B3CAE8","#2077B5","#FF7F0E","#FFBD7E","#37A537","#98DF8B","#D52728","#FF9E9C","#9467BC","#C5B1D5","#905B52","#C59B96","#E276C1","#F7B9D4","#7E7E7E","#C6C6C6","#C5C53D","#CACA52","#DBDB8D","#939B6F"))
cds.ag.differentiation <- aggregate(.~seurat_clusters,pt[,c("seurat_clusters","Pseudotime")],mean)
cds.ag.differentiation <- cds.ag.differentiation[order(cds.ag.differentiation$Pseudotime,decreasing = F),]
pt$seurat_clusters <- factor(pt$seurat_clusters,levels = cds.ag.differentiation$seurat_clusters)

p <- ggplot(pt,aes(seurat_clusters,divergence))+geom_boxplot(aes(fill=seurat_clusters),alpha=0.5,outlier.shape = NA)+
  geom_jitter(position=position_jitter(width=0.3,height=0.2),aes(colour=seurat_clusters),alpha=0.9,size=0.3)+
  scale_color_manual(values = df.color[as.character(cds.ag.differentiation$seurat_clusters),"color"])+
  scale_fill_manual(values = df.color[as.character(cds.ag.differentiation$seurat_clusters),"color"])+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",title = "",x="Cell clusters",y="AB divergence",fill="")+
  guides(color='none',fill='none')
ggsave("ABD_divergence_in_allclusters.pdf",p,width = 6.8,height = 4.6,dpi = 300)