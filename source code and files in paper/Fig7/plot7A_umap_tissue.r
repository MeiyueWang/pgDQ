library(openxlsx)
library(dplyr)
library(ggplot2)
library(pheatmap)
dinfo <- read.xlsx("info.xlsx",sheet = 1)
rownames(dinfo) <- dinfo$run_accession
pbmc <- readRDS("seurat_bulk_tissue_development.rds")
df.pbmc <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
df.pbmc$cluster <- as.factor(pbmc@meta.data$seurat_clusters)
df.pbmc$sample <- rownames(df.pbmc)
df.pbmc$tissue <- dinfo[df.pbmc$sample,"Tissue"]
df.pbmc$tissue.h <- dinfo[df.pbmc$sample,"High.level.tissue"]
df.pbmc$age <- dinfo[df.pbmc$sample,"Age"]
df.pbmc$age.h <- dinfo[df.pbmc$sample,"High.level.age"]
df.pbmc$no <- dinfo[df.pbmc$sample,"No"]
df.pbmc$no <- factor(df.pbmc$no,levels = unique(sort(df.pbmc$no)))
t1 <- colorRampPalette(c("#2077B5","white"))(length(unique(sort(df.pbmc[df.pbmc$tissue.h=="roots","no"])))+5)
t1 <- t1[c(1:(length(t1)-5))]
t2 <- colorRampPalette(c("#FF7F0E","white"))(length(unique(sort(df.pbmc[df.pbmc$tissue.h=="leaves/shoots","no"])))+5)
t2 <- t2[c(1:(length(t2)-5))]
t3 <- colorRampPalette(c("#37A537","white"))(length(unique(sort(df.pbmc[df.pbmc$tissue.h=="spike","no"])))+5)
t3 <- t3[c(1:(length(t3)-5))]
t4 <- colorRampPalette(c("#D52728","white"))(length(unique(sort(df.pbmc[df.pbmc$tissue.h=="grain","no"])))+5)
t4 <- t4[c(1:(length(t4)-5))]
dft <- data.frame(no=unique(sort(df.pbmc$no)),color=c(t1,t2,t3,t4),stringsAsFactors = F)
df.pbmc.j <- left_join(df.pbmc,dft,by=c("no"="no"))
p <- ggplot(df.pbmc,aes(UMAP_1,UMAP_2,color=no))+geom_point()+
  scale_color_manual(values = dft$color)+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "",x="UMAP 1",y="UMAP 2")+
  guides(color="none")
ggsave("tissue_UMAP_color_by_tissue.pdf",p,height = 4,width = 4,units = "in",dpi = 300)
tmpmat <- rnorm(100)
tmpmat <- matrix(tmpmat,nrow = 50,byrow = T)
pheatmap(tmpmat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#2077B5","white"))(50)) #root
pheatmap(tmpmat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#FF7F0E","white"))(50)) #leaf
pheatmap(tmpmat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#37A537","white"))(50)) #spike
pheatmap(tmpmat,cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("#D52728","white"))(50)) #grain
write.table(df.pbmc,"df.pbmc_bulk_tissue_development.xls",col.names = T,row.names = F,quote = F,sep = "\t")