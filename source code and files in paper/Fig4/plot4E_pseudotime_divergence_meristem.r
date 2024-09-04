library(ggplot2)
pt <- read.table("pseudotime_divergence_in_cells.xls",sep = "\t",stringsAsFactors = F,header = T)
p <- ggplot(pt[pt$seurat_clusters==7,],aes(Pseudotime,divergence))+geom_point(color="#2077B5")+
  geom_smooth(method = "loess",formula = y ~ x,color="#F4B183",fill="#F4B183")+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Pseudotime",y="Subgenome divergence")
ggsave("divergence_vs_pseudotime_cluster7.pdf",p,width = 5.2,height = 4.4,dpi = 300)