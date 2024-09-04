library(Seurat)
library(rdist)
library(SeuratWrappers)
library(ggplot2)

orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
CS_scRNA <- readRDS("CS_scRNA_mergereps.rds")
CS_scRNA.sel <- CS_scRNA
df <- as.data.frame(CS_scRNA.sel@assays$RNA@data)
df_A <- df[orth$cs_A,];colnames(df_A) <- paste(colnames(df_A),"A",sep = "_")
df_B <- df[orth$cs_B,];colnames(df_B) <- paste(colnames(df_B),"B",sep = "_")
df_D <- df[orth$cs_D,];colnames(df_D) <- paste(colnames(df_D),"D",sep = "_")
df_ABD <- cbind(df_A,df_B,df_D)
rownames(df_ABD) <- orth$groups
df_ABD[is.na(df_ABD)] <- 0
df.sel <- data.frame(CS_scRNA.sel$seurat_clusters,stringsAsFactors = F)
df.sel$CS_scRNA.sel.seurat_clusters <- as.character(df.sel$CS_scRNA.sel.seurat_clusters)
df.sel$cell <- gsub("-",".",rownames(df.sel))
tmpdf <- as.data.frame(matrix(c(1,0,0,0,0.5,0.5,0,1,0,0.5,0,0.5,0,0,1,0.5,0.5,0,0.33,0.33,0.33),byrow = T,ncol = 3))
f_class <- function(x){
  x <- as.numeric(x)
  tmpmat <- rbind(x[c(4:6)],tmpdf)
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

fratio <- function(x){
  c <- as.character(x["cell"])
  tmpdf2 <- df_ABD[,grep(c,colnames(df_ABD))]
  colnames(tmpdf2) <- gsub("_","",gsub(c,"",colnames(tmpdf2)))
  tmpdf2 <- tmpdf2 + 0.1
  tmpdf2$pA <- tmpdf2$A/(tmpdf2$A+tmpdf2$B+tmpdf2$D)
  tmpdf2$pB <- tmpdf2$B/(tmpdf2$A+tmpdf2$B+tmpdf2$D)
  tmpdf2$pD <- tmpdf2$D/(tmpdf2$A+tmpdf2$B+tmpdf2$D)
  tmpdf2$class <- apply(tmpdf2,1,f_class)
  #r <- nrow(tmpdf2[tmpdf2$class!="7Balanced",])/nrow(tmpdf2)
  r <- paste(tmpdf2$class,collapse=",")
}
df.sel$class <- apply(df.sel,1,fratio)
write.table(df.sel,paste("inr0","unbalanced_cell_class_all.xls",sep = "_"),row.names = F,col.names = T,quote = F,sep = "\t")
