library(ggplot2)
library(reshape2)
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
CS_scRNA <- readRDS("CS_mergereps.rds")
df <- as.data.frame(CS_scRNA@assays$RNA@data)
df.sel <- as.data.frame(t(df[c("TraesCS4A02G410800","TraesCS4B02G310900","TraesCS4D02G308800"),]))
#df.sel <- as.data.frame(t(df[c("TraesCS3A02G313500","TraesCS3B02G158400","TraesCS3D02G140800"),]))
df_out <- as.data.frame(CS_scRNA@reductions$umap@cell.embeddings)
df_out$A <- df.sel[rownames(df_out),1]
df_out$B <- df.sel[rownames(df_out),2]
df_out$D <- df.sel[rownames(df_out),3]
df_out.Z <- t(scale(t(df_out[,c("A","B","D")]),center = T,scale = T))
colnames(df_out.Z) <- c("AZ","BZ","DZ")
df_out2 <- cbind(df_out,df_out.Z)
df_out2$Anor <- ifelse(df_out2$A>6,6,df_out2$A)
df_out2$Bnor <- ifelse(df_out2$B>6,6,df_out2$B)
df_out2$Dnor <- ifelse(df_out2$D>6,6,df_out2$D)
df_out2_melt <- reshape2::melt(df_out2,id.vars=c("UMAP_1","UMAP_2"),measure.vars=c("A","B","D"))
df_out2_melt$value <- ifelse(df_out2_melt$value>5,5,df_out2_melt$value)
p <- ggplot(df_out2_melt,aes(UMAP_1,UMAP_2,color=value))+geom_point(size=0.1)+
  scale_colour_gradientn(colors = rev(c("#DD6434","white","#297AC2")))+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",title = "pg.A",x="UMAP 1",y="UMAP 2")+
  facet_grid(cols = vars(variable))
ggsave("TaTIP1.pdf",p,height = 4,width = 10,units = "in",dpi = 300)