library(rdist)
library(ggplot2)
library(ggpubr)
CS_scRNA <- readRDS("CS_mergereps.rds")
CS_scRNA.df <- as.data.frame(CS_scRNA@assays$RNA@data)
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
df_A <- CS_scRNA.df[orth$cs_A,];colnames(df_A) <- paste(colnames(df_A),"A",sep = "_")
df_B <- CS_scRNA.df[orth$cs_B,];colnames(df_B) <- paste(colnames(df_B),"B",sep = "_")
df_D <- CS_scRNA.df[orth$cs_D,];colnames(df_D) <- paste(colnames(df_D),"D",sep = "_")
df_ABD <- cbind(df_A,df_B,df_D)
rownames(df_ABD) <- orth$groups
df_ABD[is.na(df_ABD)] <- 0
CS_scRNA.cluster <- data.frame(CS_scRNA$seurat_clusters,stringsAsFactors = F)
CS_scRNA.cluster$CS_scRNA.seurat_clusters <- as.character(CS_scRNA.cluster$CS_scRNA.seurat_clusters)
CS_scRNA.cluster$cell <- rownames(CS_scRNA.cluster)
CS_scRNA.cluster$distanceAB <- apply(CS_scRNA.cluster,1,function(x){c <- as.character(x["cell"]);r <- rdist(rbind(c(df_ABD[[paste(c,"A",sep = "_")]]),c(df_ABD[[paste(c,"B",sep = "_")]])));return(r)})
CS_scRNA.cluster$distanceAD <- apply(CS_scRNA.cluster,1,function(x){c <- as.character(x["cell"]);r <- rdist(rbind(c(df_ABD[[paste(c,"A",sep = "_")]]),c(df_ABD[[paste(c,"D",sep = "_")]])));return(r)})
CS_scRNA.cluster$distanceBD <- apply(CS_scRNA.cluster,1,function(x){c <- as.character(x["cell"]);r <- rdist(rbind(c(df_ABD[[paste(c,"B",sep = "_")]]),c(df_ABD[[paste(c,"D",sep = "_")]])));return(r)})
df <- as.data.frame(CS_scRNA@reductions$umap@cell.embeddings)
df$distanceAB <- CS_scRNA.cluster[rownames(df),"distanceAB"]
df$distanceAD <- CS_scRNA.cluster[rownames(df),"distanceAD"]
df$distanceBD <- CS_scRNA.cluster[rownames(df),"distanceBD"]
df$distanceABD <- df$distanceAB+df$distanceAD+df$distanceBD
df$cluster <- CS_scRNA.cluster[rownames(df),"CS_scRNA.seurat_clusters"]
df$cluster <- factor(df$cluster,levels = as.character(c(0:19)))
p1 <- ggplot(df,aes(distanceABD,distanceAB))+geom_point(size=0.01)+
  geom_smooth(method = "lm",formula = y~x)+
  xlim(15,150)+ylim(15,150)+
  geom_abline(slope = 1)+
  #scale_colour_gradientn(colors = rev(c("#E96479","white","#297AC2","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="ABD divergence",y="AB divergence",color="",title = "")+
  guides(color=F)
p2 <- ggplot(df,aes(distanceABD,distanceAD))+geom_point(size=0.01)+
  geom_smooth(method = "lm",formula = y~x)+
  xlim(15,150)+ylim(15,150)+
  geom_abline(slope = 1)+
  #scale_colour_gradientn(colors = rev(c("#E96479","white","#297AC2","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="ABD divergence",y="AD divergence",color="",title = "")+
  guides(color=F)
p3 <- ggplot(df,aes(distanceABD,distanceBD))+geom_point(size=0.01)+
  geom_smooth(method = "lm",formula = y~x)+
  xlim(15,150)+ylim(15,150)+
  geom_abline(slope = 1)+
  #scale_colour_gradientn(colors = rev(c("#E96479","white","#297AC2","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="ABD divergence",y="BD divergence",color="",title = "")+
  guides(color=F)
ggsave("ABD_div_vs_AB.pdf",p1,width = 3.2,height = 2.9,dpi = 300)
ggsave("ABD_div_vs_AD.pdf",p2,width = 3.2,height = 2.9,dpi = 300)
ggsave("ABD_div_vs_BD.pdf",p3,width = 3.2,height = 2.9,dpi = 300)