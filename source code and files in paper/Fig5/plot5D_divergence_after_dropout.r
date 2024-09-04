library(Seurat)
library(rdist)
library(SeuratWrappers)
library(ggplot2)
dc <- read.table("divGene_out.20240412.xls",sep = "\t",stringsAsFactors = F)
dc <- dc[order(dc$V3,decreasing = T),]
dc$no <- c(1:nrow(dc))
dc$nor <- log(2/(1+dc$V3))
dc.sel <- dc[c(1:round(nrow(dc)*0.025)),"V1"]
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
dc$A <- orth[dc$V1,"cs_A"]
dc$B <- orth[dc$V1,"cs_B"]
dc$D <- orth[dc$V1,"cs_D"]
CS_scRNA <- readRDS("CS_mergereps.rds")
df <- as.data.frame(CS_scRNA@assays$RNA@data)
df_A <- df[orth$cs_A,];colnames(df_A) <- paste(colnames(df_A),"A",sep = "_")
df_B <- df[orth$cs_B,];colnames(df_B) <- paste(colnames(df_B),"B",sep = "_")
df_D <- df[orth$cs_D,];colnames(df_D) <- paste(colnames(df_D),"D",sep = "_")
df_ABD <- cbind(df_A,df_B,df_D)
rownames(df_ABD) <- orth$groups
df_ABD[is.na(df_ABD)] <- 0
df.sel <- data.frame(CS_scRNA$seurat_clusters,stringsAsFactors = F)
df.sel$CS_scRNA.seurat_clusters <- as.character(df.sel$CS_scRNA.seurat_clusters)
df.sel$cell <- rownames(df.sel)
df.sel$sc_distance0 <- apply(df.sel,1,function(x){c <- as.character(x["cell"]);r <- rdist(rbind(c(df_ABD[[paste(c,"A",sep = "_")]]),c(df_ABD[[paste(c,"B",sep = "_")]])))+rdist(rbind(c(df_ABD[[paste(c,"A",sep = "_")]]),c(df_ABD[[paste(c,"D",sep = "_")]])))+rdist(rbind(c(df_ABD[[paste(c,"B",sep = "_")]]),c(df_ABD[[paste(c,"D",sep = "_")]])));return(r)})
df_ABD.ko <- df_ABD
df_ABD.ko[dc.sel,] <- 0
df.sel$sc_distanceko <- apply(df.sel,1,function(x){c <- as.character(x["cell"]);r <- rdist(rbind(c(df_ABD.ko[[paste(c,"A",sep = "_")]]),c(df_ABD.ko[[paste(c,"B",sep = "_")]])))+rdist(rbind(c(df_ABD.ko[[paste(c,"A",sep = "_")]]),c(df_ABD.ko[[paste(c,"D",sep = "_")]])))+rdist(rbind(c(df_ABD.ko[[paste(c,"B",sep = "_")]]),c(df_ABD.ko[[paste(c,"D",sep = "_")]])));return(r)})
df_out <- as.data.frame(CS_scRNA@reductions$umap@cell.embeddings)
df_out$distanceWT <- df.sel[rownames(df_out),"sc_distance0"]
df_out$distanceKO <- df.sel[rownames(df_out),"sc_distanceko"]
df_out_melt <- reshape2::melt(df_out,id.vars=c("UMAP_1","UMAP_2"),measure.vars=c("distanceWT","distanceKO"))
p <- ggplot(df_out_melt,aes(UMAP_1,UMAP_2,color=value))+geom_point(size=0.1)+
  scale_colour_gradientn(colors = rev(c("#DD6434","white","#297AC2","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "WT",x="UMAP 1",y="UMAP 2")+
  facet_grid(cols = vars(variable))
ggsave("umap_of_divergence_after_drop_out.pdf",p,width = 6,height = 3,dpi = 300)