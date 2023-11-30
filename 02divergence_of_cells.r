library(Seurat)
library(rdist)
args <- commandArgs(T)
scRNAin <- args[1]
orthin <- args[2]
orth <- read.table(orthin,sep = "\t",stringsAsFactors = F,header = T)
rownames(orth) <- orth[[1]]
scRNA <- readRDS(scRNAin)
df <- as.data.frame(scRNA@assays$RNA@data)
df_allsubG <- c()
for (i in colnames(orth)[2:ncol(orth)]) {
  tmpdf <- df[orth[[i]],]
  tmpdf[is.na(tmpdf)] <- 0
  rownames(tmpdf) <- orth[[1]]
  df_allsubG <- rbind(df_allsubG,tmpdf)
}
df_allsubG <- as.data.frame(df_allsubG)
list_allsubG <- as.list(df_allsubG)
list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = nrow(orth),byrow = F)))))})
list_allsubG_r.df <- data.frame(cells=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
write.table(list_allsubG_r.df,"pgDQ_subgenome_divergence_of_all_cells.xls",sep = "\t",col.names = F,row.names = F,quote = F)