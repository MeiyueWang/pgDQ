library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
scRNA_uniq <- readRDS(paste("CS_mergereps.rds",sep = ""))
scRNA_EM <- readRDS(paste("CS_mergereps","EM",".rds",sep = ""))
scRNA_PropUnique <- readRDS(paste("CS_mergereps","PropUnique",".rds",sep = ""))
scRNA_Rescue <- readRDS(paste("CS_mergereps","Rescue",".rds",sep = ""))
scRNA_Uniform <- readRDS(paste("CS_mergereps","Uniform",".rds",sep = ""))
df_cell <- unique(sort(c(gsub("-1","",colnames(scRNA_uniq)),
                         colnames(scRNA_EM),
                         colnames(scRNA_PropUnique),
                         colnames(scRNA_Rescue),
                         colnames(scRNA_Uniform))))
df_cell <- data.frame(df_cell,stringsAsFactors = F)
for (i in c("uniq","EM","PropUnique","Rescue","Uniform")) {
  CS_scRNA <- get(paste("scRNA",i,sep = "_"))
  CS_scRNA.df <- as.data.frame(CS_scRNA@assays$RNA@data)
  orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
  orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
  orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
  orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
  rownames(orth) <- orth$groups
  dfr.A <- CS_scRNA.df[orth$cs_A,];dfr.A[is.na(dfr.A)] <- 0;rownames(dfr.A) <- orth$groups
  dfr.B <- CS_scRNA.df[orth$cs_B,];dfr.B[is.na(dfr.B)] <- 0;rownames(dfr.B) <- orth$groups
  dfr.D <- CS_scRNA.df[orth$cs_D,];dfr.D[is.na(dfr.D)] <- 0;rownames(dfr.D) <- orth$groups
  df_allsubG <- rbind(dfr.A,dfr.B,dfr.D)
  df_allsubG <- as.data.frame(df_allsubG)
  list_allsubG <- as.list(df_allsubG)
  list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = nrow(orth),byrow = F)))))})
  list_allsubG_r.df <- data.frame(samples=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
  list_allsubG_r.df$cell <- gsub("-1","",list_allsubG_r.df$samples)
  df_cell <- left_join(df_cell,list_allsubG_r.df,by=c("df_cell"="cell"))
}
rownames(df_cell) <- df_cell$df_cell
df_cell <- df_cell[,c("divergence.x","divergence.y","divergence.x.x","divergence.y.y","divergence")]
colnames(df_cell) <- c("cellranger","STARsolo_EM","STARsolo_PropUnique","STARsolo_Rescue","STARsolo_Uniform")
df_cell$cells <- rownames(df_cell)
#write.table(df_cell,"divergence_of_all_mapping_algorithm.xls",quote = F,sep = "\t",row.names = F,col.names = T)
p1 <- ggplot(na.omit(df_cell[,c("cellranger","STARsolo_EM")]),aes(cellranger,STARsolo_EM))+geom_point(size=0.5)+
  geom_smooth(method = "lm",formula = y~x)+
  stat_cor(method = "pearson",color="red",size=5)+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Cell Ranger",y="STARsolo EM")
p2 <- ggplot(na.omit(df_cell[,c("cellranger","STARsolo_Uniform")]),aes(cellranger,STARsolo_Uniform))+geom_point(size=0.5)+
  geom_smooth(method = "lm",formula = y~x)+
  stat_cor(method = "pearson",color="red",size=5)+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Cell Ranger",y="STARsolo Uniform")
p3 <- ggplot(na.omit(df_cell[,c("cellranger","STARsolo_PropUnique")]),aes(cellranger,STARsolo_PropUnique))+geom_point(size=0.5)+
  geom_smooth(method = "lm",formula = y~x)+
  stat_cor(method = "pearson",color="red",size=5)+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Cell Ranger",y="STARsolo PropUnique")
p4 <- ggplot(na.omit(df_cell[,c("cellranger","STARsolo_Rescue")]),aes(cellranger,STARsolo_Rescue))+geom_point(size=0.5)+
  geom_smooth(method = "lm",formula = y~x)+
  stat_cor(method = "pearson",color="red",size=5)+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Cell Ranger",y="STARsolo Rescue")


ggsave("cellranger_vs_EM.pdf",p1,width = 3.2,height = 2.9,dpi = 300)
ggsave("cellranger_vs_Uniform.pdf",p2,width = 3.2,height = 2.9,dpi = 300)
ggsave("cellranger_vs_PropUnique.pdf",p3,width = 3.2,height = 2.9,dpi = 300)
ggsave("cellranger_vs_Rescue.pdf",p4,width = 3.2,height = 2.9,dpi = 300)


t1 <- na.omit(df_cell[,c("cellranger","STARsolo_EM")])
t2 <- na.omit(df_cell[,c("cellranger","STARsolo_Uniform")])
t3 <- na.omit(df_cell[,c("cellranger","STARsolo_PropUnique")])
t4 <- na.omit(df_cell[,c("cellranger","STARsolo_Rescue")])
cor(t1$cellranger,t1$STARsolo_EM)
cor(t2$cellranger,t2$STARsolo_Uniform)
cor(t3$cellranger,t3$STARsolo_PropUnique)
cor(t4$cellranger,t4$STARsolo_Rescue)