library(Seurat)
library(rdist)
library(ggplot2)
library(ggtern)
tmpdf <- as.data.frame(matrix(c(1,0,0,0,0.5,0.5,0,1,0,0.5,0,0.5,0,0,1,0.5,0.5,0,0.33,0.33,0.33),byrow = T,ncol = 3))
f_class <- function(x){
  x <- as.numeric(x)
  tmpmat <- rbind(x[c(1:3)],tmpdf)
  dist_df <- as.data.frame(as.matrix(rdist(tmpmat)))[-1,]
  minpos <- which.min(dist_df[,1])
  if(minpos==1){
    c <- "1A dominant"
  }else if(minpos==2){
    c <- "2A suppressed"
  }else if(minpos==3){
    c <- "3B dominant"
  }else if(minpos==4){
    c <- "4B suppressed"
  }else if(minpos==5){
    c <- "5D dominant"
  }else if(minpos==6){
    c <- "6D suppressed"
  }else{
    c <- "7Balanced"
  }
  return(c)
}
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
root_inh <- read.table("CS_root_expression.xls",stringsAsFactors = F,sep = "\t",row.names = 1,header = T)
root_inh_1.A <- root_inh[orth$cs_A,"CS_root_rep1"]+0.00001
root_inh_1.B <- root_inh[orth$cs_B,"CS_root_rep1"]+0.00001
root_inh_1.D <- root_inh[orth$cs_D,"CS_root_rep1"]+0.00001
root_inh_1.ABD <- data.frame(root_inh_1.A,root_inh_1.B,root_inh_1.D,stringsAsFactors = F)
rownames(root_inh_1.ABD) <- orth$groups
root_inh_1.ABD$Apercent <- root_inh_1.ABD$root_inh_1.A/(root_inh_1.ABD$root_inh_1.A+root_inh_1.ABD$root_inh_1.B+root_inh_1.ABD$root_inh_1.D)
root_inh_1.ABD$Bpercent <- root_inh_1.ABD$root_inh_1.B/(root_inh_1.ABD$root_inh_1.A+root_inh_1.ABD$root_inh_1.B+root_inh_1.ABD$root_inh_1.D)
root_inh_1.ABD$Dpercent <- root_inh_1.ABD$root_inh_1.D/(root_inh_1.ABD$root_inh_1.A+root_inh_1.ABD$root_inh_1.B+root_inh_1.ABD$root_inh_1.D)
root_inh_1.ABD$class <- apply(root_inh_1.ABD[,c("Apercent","Bpercent","Dpercent")],1,f_class)
root_inh_1.ABD <- root_inh_1.ABD[order(root_inh_1.ABD$class),]
p <- ggplot()+coord_tern()+geom_point(data=subset(root_inh_1.ABD,class=="1Balanced"),aes(Apercent,Bpercent,Dpercent),color="#FDF7B8",size=0.2)+
  geom_point(data=subset(root_inh_1.ABD,class!="1Balanced"),aes(Apercent,Bpercent,Dpercent,color=class),size=0.2)+
  theme_bw(base_size = 20)+
  scale_color_manual(values = c("#5BA346","#A1D298","#4F2E7F","#30BBFF","#F9962D","#FFCB22","#9DA1A0"),labels=c("A dominant","A suppressed","B dominant","B suppressed","D dominant","D suppressed","Balanced"))+
  labs(x="A",y="B",z="D",color="")+
  theme(axis.text = element_text(colour = "black",size = 20),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 2,colour = "black"),
        axis.title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0,0,0,0),"cm"))
  guides(color="none")
ggsave("triangle_plot_bulk_root.pdf",p,width = 3.3,height = 3.3,dpi = 300)