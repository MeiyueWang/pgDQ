library(rdist)
library(ggplot2)
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

mono <- read.table("pgDQ_monocle_pseudotime.xls",header = T)
rownames(mono) <- mono$cells
df$Pseudotime <- mono[rownames(df),"Pseudotime"]

df.color <- data.frame(no=c(0:19),color=c("#B3CAE8","#2077B5","#FF7F0E","#FFBD7E","#37A537","#98DF8B","#D52728","#FF9E9C","#9467BC","#C5B1D5","#905B52","#C59B96","#E276C1","#F7B9D4","#7E7E7E","#C6C6C6","#C5C53D","#CACA52","#DBDB8D","#939B6F"))
rownames(df.color) <- df.color$no

cds.ag.differentiation <- aggregate(.~cluster,df[,c("cluster","Pseudotime")],mean)
cds.ag.differentiation <- cds.ag.differentiation[order(cds.ag.differentiation$Pseudotime,decreasing = F),]
df$cluster <- factor(df$cluster,levels = cds.ag.differentiation$cluster)

p1 <- ggplot(df,aes(cluster,distanceAB))+geom_boxplot(aes(fill=cluster),alpha=0.5,outlier.shape = NA)+
  geom_jitter(position=position_jitter(width=0.3,height=0.2),aes(colour=cluster),alpha=0.9,size=0.3)+
  scale_color_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  scale_fill_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",title = "",x="Cell clusters",y="AB divergence",fill="")+
  guides(color='none',fill='none')
p2 <- ggplot(df,aes(cluster,distanceAD))+geom_boxplot(aes(fill=cluster),alpha=0.5,outlier.shape = NA)+
  geom_jitter(position=position_jitter(width=0.3,height=0.2),aes(colour=cluster),alpha=0.9,size=0.3)+
  scale_color_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  scale_fill_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",title = "",x="Cell clusters",y="AD divergence",fill="")+
  guides(color='none',fill='none')
p3 <- ggplot(df,aes(cluster,distanceBD))+geom_boxplot(aes(fill=cluster),alpha=0.5,outlier.shape = NA)+
  geom_jitter(position=position_jitter(width=0.3,height=0.2),aes(colour=cluster),alpha=0.9,size=0.3)+
  scale_color_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  scale_fill_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",title = "",x="Cell clusters",y="BD divergence",fill="")+
  guides(color='none',fill='none')
p4 <- ggplot(df,aes(cluster,distanceABD))+geom_boxplot(aes(fill=cluster),alpha=0.5,outlier.shape = NA)+
  geom_jitter(position=position_jitter(width=0.3,height=0.2),aes(colour=cluster),alpha=0.9,size=0.3)+
  scale_color_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  scale_fill_manual(values = df.color[as.character(cds.ag.differentiation$cluster),"color"])+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",title = "",x="Cell clusters",y="ABD divergence",fill="")+
  guides(color='none',fill='none')
ggsave("AB_divergence_in_allclusters.pdf",p1,width = 6.8,height = 4.6,dpi = 300)
ggsave("AD_divergence_in_allclusters.pdf",p2,width = 6.8,height = 4.6,dpi = 300)
ggsave("BD_divergence_in_allclusters.pdf",p3,width = 6.8,height = 4.6,dpi = 300)
ggsave("ABD_divergence_in_allclusters.pdf",p4,width = 6.8,height = 4.6,dpi = 300)


ABc <- cor(df$distanceAB,df$pt,method = "spearman")
ADc <- cor(df$distanceAD,df$pt,method = "spearman")
BDc <- cor(df$distanceBD,df$pt,method = "spearman")
ABDc <- cor(df$distanceABD,df$pt,method = "spearman")
df2 <- data.frame(sample=c("AB","AD","BD","ABD"),cor=c(ABc,ADc,BDc,ABDc))
df2$sample <- factor(df2$sample,levels = c("AB","AD","BD","ABD"))
ggplot(df2,aes(sample,cor))+geom_bar(stat = "identity",fill="steelblue")+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Comparison",y="Correlation between \ndistance and pseudotime")