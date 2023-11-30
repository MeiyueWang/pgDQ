library(Seurat)
library(monocle)

args <- commandArgs(T)
scRNAin <- args[1]

scRNA <- readRDS(scRNAin)
pbmc <- scRNA
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
save(cds,file="pgDQ_monocle.RData")
p1 <- plot_cell_trajectory(cds, color_by = "seurat_clusters")
ggsave("pgDQ_monocle_trajectory_colored_by_cell_clusters.png",p1,width = 5,height = 5,dpi = 300)
cds <- orderCells(cds)
p2 <- plot_cell_trajectory(cds)
ggsave("pgDQ_monocle_trajectory_colored_by_state.png",p2,width = 5,height = 5,dpi = 300)