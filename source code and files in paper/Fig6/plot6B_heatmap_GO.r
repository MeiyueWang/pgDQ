library(ggplot2)
library(pheatmap)
library(data.table)
library(dplyr)
library(rdist)
data <- read.table("scRNA_GO_150_1000.Gsea.1712991672305/gsea_report_for_IS_1712991672305.tsv",sep = "\t",stringsAsFactors = F,header = T)
GOinfo <- fread("GO/GO.info",sep = "\t",stringsAsFactors = F,header = F)
data.j <- left_join(data,GOinfo,by=c("GS.br..follow.link.to.MSigDB"="V1"))

orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
GO <- read.table("GO/GO.anno",sep = "\t",stringsAsFactors = F)
cds.data <- read.table("pgDQ_monocle_pseudotime.xls",sep = "\t",header = T,stringsAsFactors = F)
CS_scRNA <- readRDS("CS_mergereps.rds")
df <- as.data.frame(CS_scRNA@assays$RNA@data)
df_A <- df[orth$cs_A,];rownames(df_A) <- orth$groups;df_A[is.na(df_A)] <- 0
df_B <- df[orth$cs_B,];rownames(df_B) <- orth$groups;df_B[is.na(df_B)] <- 0
df_D <- df[orth$cs_D,];rownames(df_D) <- orth$groups;df_D[is.na(df_D)] <- 0
gsea_IS <- read.table("scRNA_GO_150_1000.Gsea.1712991672305/gsea_report_for_IS_1712991672305.tsv",sep = "\t",stringsAsFactors = F,header = T)
gsea_IS <- gsea_IS[gsea_IS$NES>=1,]
gsea_zero <- read.table("scRNA_GO_150_1000.Gsea.1712991672305/gsea_report_for_zero_1712991672305.tsv",sep = "\t",stringsAsFactors = F,header = T)
gsea_zero <- gsea_zero[gsea_zero$NES<= -1,]
gsea_IS_GO <- c()
for (i in c(gsea_IS$NAME[c(1,2,6,7,8,9,10,12,15,16)])) {
  GO_tG <- GO[grep(i,GO$V2,ignore.case = T),"V1"]
  GO_tG_og <- unique(sort(c(orth[orth$cs_A%in%GO_tG,"groups"],orth[orth$cs_B%in%GO_tG,"groups"],orth[orth$cs_D%in%GO_tG,"groups"])))
  gsea_IS_GO <- c(gsea_IS_GO,i)
  df_A.sel <- df_A[GO_tG_og,]
  df_B.sel <- df_B[GO_tG_og,]
  df_D.sel <- df_D[GO_tG_og,]
  df_allsubG <- rbind(df_A.sel,df_B.sel,df_D.sel)
  df_allsubG <- as.data.frame(df_allsubG)
  list_allsubG <- as.list(df_allsubG)
  list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = length(GO_tG_og),byrow = F)))))})
  list_allsubG_r.df <- data.frame(samples=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
  rownames(list_allsubG_r.df) <- list_allsubG_r.df$samples
  cds.data[[paste(i,"div",sep = "_")]] <- list_allsubG_r.df[cds.data$cells,"divergence"]
}
gsea_zero_GO <-  c()
for (i in c(gsea_zero$NAME[c(1:10)])) {
  GO_tG <- GO[grep(i,GO$V2,ignore.case = T),"V1"]
  GO_tG_og <- unique(sort(c(orth[orth$cs_A%in%GO_tG,"groups"],orth[orth$cs_B%in%GO_tG,"groups"],orth[orth$cs_D%in%GO_tG,"groups"])))
  gsea_zero_GO <- c(gsea_zero_GO,i)
  df_A.sel <- df_A[GO_tG_og,]
  df_B.sel <- df_B[GO_tG_og,]
  df_D.sel <- df_D[GO_tG_og,]
  df_allsubG <- rbind(df_A.sel,df_B.sel,df_D.sel)
  df_allsubG <- as.data.frame(df_allsubG)
  list_allsubG <- as.list(df_allsubG)
  list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = length(GO_tG_og),byrow = F)))))})
  list_allsubG_r.df <- data.frame(samples=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
  rownames(list_allsubG_r.df) <- list_allsubG_r.df$samples
  cds.data[[paste(i,"div",sep = "_")]] <- list_allsubG_r.df[cds.data$cells,"divergence"]
}

mat.sel <- cds.data[,c("Pseudotime",paste(c(gsea_IS_GO,gsea_zero_GO),"div",sep = "_"))]
mat.sel <- mat.sel[order(mat.sel$Pseudotime),]
mat.sel <- mat.sel[,grep("div",colnames(mat.sel))]
mat.sel$sum <- apply(mat.sel,1,sum)
mat.sel.nor <- apply(mat.sel[mat.sel$sum>0,-ncol(mat.sel)],2,function(x){scale(x,center = T,scale = T)})
mat.sel.nor <- apply(mat.sel.nor,2,function(x){ifelse(x>1.5,1.5,ifelse(x< -1.5,-1.5,x))})
mat.sel.nor <- as.data.frame(mat.sel.nor)


pdf("GO_pathway_triads_divergence.pdf",width = 5,height = 2)
pheatmap(t(mat.sel.nor[c(1:4935),paste(gsea_IS_GO,"div",sep = "_")]),cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,color = colorRampPalette(c("#FFFFCC","#FFF1AA","#FEE187","#FEC965","#FEAB49","#ED2F22","#B00026"))(50),border_color = NA)
dev.off()
#### "#5E50A1","#A6CE55","#FDF7B8","#F15B2B","#A51E25"
####"black","#FAC402"
####c("#FFFFCC","#FFF1AA","#FEE187","#FEC965","#FEAB49","#ED2F22","#B00026")
pdf("GO_pathway_triads_divergence_bg.pdf",width = 5,height = 1.5)
pheatmap(t(mat.sel.nor[c(1:4935),paste(gsea_zero_GO,"div",sep = "_")]),cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,color = colorRampPalette(c("#297AC2","white","#E96479"))(50),border_color = NA)
dev.off()

GOinfo <- fread("GO/GO.info",sep = "\t",stringsAsFactors = F,header = F)
gsea_IS_GO.df <- data.frame(GO=gsea_IS_GO,no=c(1:length(gsea_IS_GO)),stringsAsFactors = F)
gsea_IS_GO.df.j <- left_join(gsea_IS_GO.df,GOinfo,by=c("GO"="V1"))
gsea_zero_GO.df <- data.frame(GO=gsea_zero_GO,no=c(1:length(gsea_zero_GO)),stringsAsFactors = F)
gsea_zero_GO.df.j <- left_join(gsea_zero_GO.df,GOinfo,by=c("GO"="V1"))