library(monocle)
pbmc <- readRDS("seurat_bulk_tissue_development.rds")
umi <- as(as.matrix(pbmc@assays$RNA@counts), 'sparseMatrix')
pData <- pbmc@meta.data
pData$celltype <- pbmc@active.ident
fData <- data.frame(
  gene_short_name = row.names(pbmc),
  row.names = row.names(pbmc)
)
pd <- new('AnnotatedDataFrame', data=pData)
fd <- new('AnnotatedDataFrame', data=fData)
cds <- newCellDataSet(
  umi,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.1,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
deg.cluster <- FindAllMarkers(pbmc)
express_genes <- subset(deg.cluster, p_val_adj < 0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
plot_cell_trajectory(cds, color_by = "seurat_clusters")
cds <- orderCells(cds,root_state = 5)
plot_cell_trajectory(cds,color_by = "Pseudotime")
cds.data <- pData(cds)
save(cds,file="monocle_bulk_tissue_development.RData")
write.table(cds.data, "monocle_bulk_tissue_development_cds.data.txt",sep="\t",col.names=T,row.names=F,quote=F )