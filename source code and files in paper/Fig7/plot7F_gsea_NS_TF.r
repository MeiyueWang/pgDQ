library(ggplot2)
library(cowplot)
library(data.table)
library(dplyr)
df_IS <- read.table("df_IS_tissue.txt",sep = "\t",stringsAsFactors = F,header = T)
df_IS$no <- df_IS$no <- c(1:nrow(df_IS))
## Load .edb data
gsea.edb <- read.delim("tissue_TF_20231217.Gsea.1702819019427/edb/results.edb",header = FALSE,stringsAsFactors = FALSE)
gsea.edb <- unlist(gsea.edb)

fsel <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  NES <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NES="))[2]," NP="))[1])
  NP <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NP="))[2]," FDR="))[1])
  FDR <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"FDR="))[2]," FWER="))[1])
  RND_ES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"RND_ES="))[2]," HIT_INDICES="))[1]," ")))
  HIT_INDICES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"HIT_INDICES="))[2]," ES_PROFILE="))[1]," ")))
  ES_PROFILE <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"ES_PROFILE="))[2]," RANK_AT_ES="))[1]," ")))
}
fNES <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  NES <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NES="))[2]," NP="))[1])
}
fNP <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  NP <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NP="))[2]," FDR="))[1])
}
fFDR <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  FDR <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"FDR="))[2]," FWER="))[1])
}
fRND_ES <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  RND_ES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"RND_ES="))[2]," HIT_INDICES="))[1]," ")))
}
fHIT_INDICES <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  HIT_INDICES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"HIT_INDICES="))[2]," ES_PROFILE="))[1]," ")))
}
fES_PROFILE <- function(gene.set){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  ES_PROFILE <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"ES_PROFILE="))[2]," RANK_AT_ES="))[1]," ")))
}
fplot <- function(gene.set,color){
  sel <- gsea.edb[grep(gene.set,gsea.edb)]
  NES <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NES="))[2]," NP="))[1])
  NP <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"NP="))[2]," FDR="))[1])
  FDR <- as.numeric(unlist(strsplit(unlist(strsplit(sel,"FDR="))[2]," FWER="))[1])
  RND_ES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"RND_ES="))[2]," HIT_INDICES="))[1]," ")))
  HIT_INDICES <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"HIT_INDICES="))[2]," ES_PROFILE="))[1]," ")))
  ES_PROFILE <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(sel,"ES_PROFILE="))[2]," RANK_AT_ES="))[1]," ")))
  df.ES <- data.frame(ES=ES_PROFILE)
  df.ES$no <- c(1:nrow(df.ES))
  p1 <- ggplot(df.ES,aes(no,ES))+geom_line(linewidth=1,color=color)+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text = element_text(colour = "black",size = 14),
          axis.ticks = element_line(size = 1,colour = "black"),
          axis.title = element_text(size = 14),
          plot.margin = margin(5.5,5.5,0,5.5))+
    labs(x="",y="Enrichment score (ES)")+
    scale_x_discrete(expand = c(0,0),position = "top")+
    guides(color="none")
  df.hit <- data.frame(hit=HIT_INDICES)
  #df.hit2 <- data.frame(c(1:49282),stringsAsFactors = F)
  #df.hit2$tmp <- ifelse(df.hit2$c.1.49282.%in%df.hit$hit,1,0)
  df.hit$tmp <- 1
  p2 <- ggplot(df.hit,aes(hit,1,fill=tmp))+geom_tile(color=color)+
    theme_bw(base_size = 14)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(colour = "black",size = 14),
          axis.ticks = element_line(size = 1,colour = "black"),
          axis.title = element_text(size = 14),
          plot.margin = margin(0, 5.5,5.5, 5.5))+
    labs(x="",y="")+
    scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
    guides(fill=F)
  p3 <- ggplot(df_IS,aes(no,IS))+geom_bar(stat = "identity",fill="grey")+
    theme_bw(base_size = 14)+
    theme(axis.text = element_text(colour = "black",size = 14),
          axis.ticks = element_line(size = 1,colour = "black"),
          axis.title = element_text(size = 14,hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(-30, 5.5, 5.5, 5.5))+
    labs(x="",y= "normalized score")+
    scale_x_discrete(expand = c(0,0))
  p <- plot_grid(p1,p2,p3,align = "v",nrow = 3,rel_heights = c(3,0.9,1.5),axis = "l")
  ggsave(paste("gsea_TF_",gsub(":","_",gene.set),"v2.pdf",sep = ""),p,width = 4,height = 4,dpi = 300)
}
fplot("MYB-2A-2","#F37258")
fplot("AP2-ERF-5D-1","#625E9B")
fplot("NAC-2A-3","#F38B97")
fplot("LOB-3B-1","#72B1DE")
fplot("B3-5B-1","#ABC0E4")
