library(Seurat)
library(monocle)
div <- read.table("divergence_bulk_tissue_development.txt",sep = "\t",stringsAsFactors = F,header = T)
rownames(div) <- div$samples
df <- read.table("all_sample_lnc_count.GeneV1.1.20200901.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
dfr <- data.frame(gene=rownames(df),stringsAsFactors = F)
i <- "Development"
dtmp <- fread(paste(i,"_tpm.tsv",sep = ""),stringsAsFactors = F,header = T,data.table = F)
dfr <- left_join(dfr,dtmp,by=c("gene"="gene"))
rownames(dfr) <- dfr$gene
dfr <- dfr[,-1]
pbmc <- CreateSeuratObject(counts = dfr,project = "SeuratProject", assay = "RNA",min.cells = 0, min.features = 0, names.field = 1,names.delim = "_", meta.data = NULL)
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize")
pbmc <- FindVariableFeatures(pbmc,selection.method = "mvp")
pbmc <- ScaleData(pbmc,features = rownames(pbmc))
pbmc <- RunPCA(pbmc,npcs = 100)
pbmc <- FindNeighbors(pbmc,dims = 1:50)
pbmc <- FindClusters(pbmc,resolution = 0.5)
pbmc <- RunTSNE(pbmc, dims = 1:50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
saveRDS(pbmc,file = "seurat_bulk_tissue_development.rds")