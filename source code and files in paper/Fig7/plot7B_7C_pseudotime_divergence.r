library(data.table)
library(ggplot2)
df.pbmc <- fread("df.pbmc_bulk_tissue_development.xls",data.table = F)
rownames(df.pbmc) <- df.pbmc$sample
div <- read.table("divergence_bulk_tissue_development.txt",sep = "\t",stringsAsFactors = F,header = T)
rownames(div) <- div$samples
df.pbmc$divergence <- div[df.pbmc$sample,"divergence"]
cds.data <- read.table("monocle_bulk_tissue_development_cds.data.txt",sep = "\t",stringsAsFactors = F,header = T)
rownames(cds.data) <- cds.data$sample
df.pbmc$pseudotime <- cds.data[df.pbmc$sample,"Pseudotime"]
p1 <- ggplot(df.pbmc,aes(UMAP_1,UMAP_2,color=pseudotime))+geom_point()+
  scale_colour_gradientn(colors = rev(c("#E96479","white","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))
p2 <- ggplot(df.pbmc,aes(UMAP_1,UMAP_2,color=divergence))+geom_point()+
  scale_colour_gradientn(colors = rev(c("#E96479","white","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))
ggsave("tissue_UMAP_color_by_pseudotime.pdf",p1,height = 4,width = 5,units = "in",dpi = 300)
ggsave("tissue_UMAP_color_by_divergence.pdf",p2,height = 4,width = 5,units = "in",dpi = 300)