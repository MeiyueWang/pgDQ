library(Seurat)
library(rdist)
library(ggtern)
CS_scRNA <- readRDS("CS_mergereps.rds")
CS_scRNA.df <- as.data.frame(CS_scRNA@assays$RNA@data)
df_A <- CS_scRNA.df[orth$cs_A,];colnames(df_A) <- paste(colnames(df_A),"A",sep = "_")
df_B <- CS_scRNA.df[orth$cs_B,];colnames(df_B) <- paste(colnames(df_B),"B",sep = "_")
df_D <- CS_scRNA.df[orth$cs_D,];colnames(df_D) <- paste(colnames(df_D),"D",sep = "_")
df_ABD <- cbind(df_A,df_B,df_D)
rownames(df_ABD) <- orth$groups
df_ABD <- na.omit(df_ABD)
onecell <- df_ABD[,grep("AAACCCACATGACAGG-1_1",colnames(df_ABD))]
onecell <- onecell + 0.00001
onecell$Apercent <- onecell$`AAACCCACATGACAGG-1_1_A`/(onecell$`AAACCCACATGACAGG-1_1_A` + onecell$`AAACCCACATGACAGG-1_1_B` + onecell$`AAACCCACATGACAGG-1_1_D`)
onecell$Bpercent <- onecell$`AAACCCACATGACAGG-1_1_B`/(onecell$`AAACCCACATGACAGG-1_1_A` + onecell$`AAACCCACATGACAGG-1_1_B` + onecell$`AAACCCACATGACAGG-1_1_D`)
onecell$Dpercent <- onecell$`AAACCCACATGACAGG-1_1_D`/(onecell$`AAACCCACATGACAGG-1_1_A` + onecell$`AAACCCACATGACAGG-1_1_B` + onecell$`AAACCCACATGACAGG-1_1_D`)
onecell$class <- apply(onecell[,c("Apercent","Bpercent","Dpercent")],1,f_class)
onecell$class <- factor(onecell$class,levels = c("1A dominant","2A suppressed","3B dominant","4B suppressed","5D dominant","6D suppressed","7Balanced"))
p <- ggplot()+coord_tern()+geom_point(data=subset(onecell,class=="1Balanced"),aes(Apercent,Bpercent,Dpercent),color="#9BC0E2",size=1)+geom_point(data=subset(onecell,class!="1Balanced"),aes(Apercent,Bpercent,Dpercent,color=class),size=1)+
  theme_bw(base_size = 30)+
  scale_color_manual(values = c("#5BA346","#A1D298","#4F2E7F","#30BBFF","#F9962D","#FFCB22","#9DA1A0"),labels=c("A dominant","A suppressed","B dominant","B suppressed","D dominant","D suppressed","Balanced"))+
  labs(x="A",y="B",z="D",color="")+
  theme(axis.text = element_text(colour = "black",size = 20),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 2,colour = "black"),
        axis.title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0,0,0,0),"cm"))+
  guides(color="none")
ggsave("triangle_plot_onecell.pdf",p,width = 3.3,height = 3.3,dpi = 300)