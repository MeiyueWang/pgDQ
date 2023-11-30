library(Seurat)

args <- commandArgs(T)
scRNAin <- args[1]
marker_genein <- args[2]

mG <- as.character(unlist(read.table(marker_genein,sep = "\t",stringsAsFactors = F,header = F)))
scRNA <- readRDS(scRNAin)
df <- as.data.frame(scRNA@reductions$umap@cell.embeddings)
df$cluster <- as.factor(scRNA@meta.data$seurat_clusters)
if (length(mG) == 1) {
  gene_df <- as.data.frame(GetAssayData(object = CS_scRNA,slot = "data")[mG,])
  colnames(gene_df) <- mG
  gene_df$exp <- gene_df[[mG]]
}else{
  gene_df <- as.data.frame(t(as.data.frame(GetAssayData(object = CS_scRNA,slot = "data")[mG,])))
  gene_df$exp <- apply(gene_df,1,sum)
}

df$exp <- gene_df[rownames(df),"exp"]
df$cluster <- as.character(df$cluster)
df.ag <- aggregate(.~cluster,df[,c("cluster","exp")],sum)
df.ag <- df.ag[order(df.ag$exp,decreasing = T),]
cluster.sel <- df.ag[1,"cluster"]
write.table(paste("initial cluster:",cluster.sel,sep = " "),"pgDQ_initial_cluster.txt",quote = F,row.names = F,col.names = F)
