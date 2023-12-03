library(Seurat)
library(rdist)
library(foreach)
library(doParallel)

args <- commandArgs(T)
scRNAin <- args[1]
orthin <- args[2]
ncores <- args[3]

cl <- makeCluster(ncores)
registerDoParallel(cl)

orth <- read.table(orthin,sep = "\t",stringsAsFactors = F,header = T)
rownames(orth) <- orth[[1]]
cds.data <- read.table("pgDQ_monocle_pseudotime.xls",sep = "\t",header = T)
rownames(cds.data) <- cds.data$cells
div <- read.table("pgDQ_subgenome_divergence_of_all_cells.xls",sep = "\t",stringsAsFactors = F,header = F)
rownames(div) <- div$V1
cds.data$divergence <- div[rownames(cds.data),"V2"]
#corr0 <- cor(cds.data$divergence,cds.data$Pseudotime,method = "spearman")
scRNA <- readRDS(scRNAin)
df <- as.data.frame(scRNA@assays$RNA@data)
t <- nrow(orth)
df_out <- foreach(j=1:t, .combine = rbind, .packages =  c('doParallel','rdist'), .errorhandling = "pass") %dopar% {
  df.ko <- df
  df.ko[as.character(orth[j,c(2:ncol(orth))]),] <- 0
  df_allsubG <- c()
  for (i in colnames(orth)[2:ncol(orth)]) {
    tmpdf <- df.ko[orth[[i]],]
    tmpdf[is.na(tmpdf)] <- 0
    rownames(tmpdf) <- orth[[1]]
    df_allsubG <- rbind(df_allsubG,tmpdf)
  }
  df_allsubG <- as.data.frame(df_allsubG)
  list_allsubG <- as.list(df_allsubG)
  list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = nrow(orth),byrow = F)))))})
  list_allsubG_r.df <- data.frame(cells=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
  list_allsubG_r.df$pt <- cds.data[list_allsubG_r.df$cells,"Pseudotime"]
  corr.t <- cor.test(list_allsubG_r.df$divergence,list_allsubG_r.df$pt,method = "spearman")
  corr <- as.numeric(corr.t$estimate)
  p <- as.numeric(corr.t$p.value)
  r <- paste(c(orth[j,1],p,corr),collapse = ";")
  return(r)
  gc()
  #write.table(r,paste(orth[j,1],"txt",sep = "."),sep = "\t",col.names = F,row.names = F,quote = F)
}
stopImplicitCluster()
df_out <- as.data.frame(df_out)
df_out2 <- as.character(unlist(strsplit(df_out$V1,split = ";")))
mat <- matrix(df_out2,byrow = T,ncol = 3)
mat <- as.data.frame(mat)
mat$V2 <- as.numeric(mat$V2)
mat$V3 <- as.numeric(mat$V3)
mat$IS <- log10(1-mat$V3)
colnames(mat) <- c("orthgroups","pvalue","corr","IS")
colnames(orth)[1] <- "orthgroups"
mat.j <- inner_join(orth,mat,by=c("orthgroups"="orthgroups"))
write.table(mat.j,"pgDQ_indispensable_scores_of_all_1v1orthologous_genes.xls",sep = "\t",col.names = T,row.names = F,quote = F)
