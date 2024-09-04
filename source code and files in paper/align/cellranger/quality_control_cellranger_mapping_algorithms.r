library(SoupX)
library(Seurat)
library(harmony)
library(DoubletFinder)
library(ggplot2)
sc_rep1.data <- load10X("CS_singlecell_rep1_cellranger_chrC_chrM")
sc_rep1.rho <- autoEstCont(sc_rep1.data)
sc_rep1.out <- adjustCounts(sc_rep1.rho)
sc_rep2.data <- load10X("CS_singlecell_rep2_cellranger_chrC_chrM")
sc_rep2.rho <- autoEstCont(sc_rep2.data)
sc_rep2.out <- adjustCounts(sc_rep2.rho)
sc_rep1 <- CreateSeuratObject(counts = sc_rep1.out,project="rep1",min.cells=3)
sc_rep2 <- CreateSeuratObject(counts = sc_rep2.out,project="rep2",min.cells=3)
sc <- merge(x=sc_rep1,y=sc_rep2)
sc[["percent.mt"]] <- PercentageFeatureSet(sc,pattern = "^Mitochondria")
#sc[["percent.cl"]] <- PercentageFeatureSet(sc,pattern = "^Chloroplast")


p1 <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
ggsave("nFeature_nCount_before_QC.pdf",p1,width = 7.46,height = 5.14,dpi = 300)
#FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5)
ggsave("nFeature_nCount_scatter_before_QC.pdf",p2,width = 4.77,height = 4.03,dpi = 300)
#FeatureScatter(sc, feature1 = "nFeature_RNA",feature2 = "percent.mt")


sc1 <- subset(sc, subset = nFeature_RNA > 500 & nFeature_RNA < 10000)
sc1 <- subset(sc1, subset = nCount_RNA > 500 & nCount_RNA < 20000)
sc1 <- subset(sc1, subset = percent.mt < 2)


p3 <- VlnPlot(sc1, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
ggsave("nFeature_nCount_after_QC.pdf",p3,width = 7.46,height = 5.14,dpi = 300)
p4 <- FeatureScatter(sc1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5)
ggsave("nFeature_nCount_scatter_after_QC.pdf",p4,width = 4.77,height = 4.03,dpi = 300)


sc.nor <- NormalizeData(sc1, normalization.method = "LogNormalize", scale.factor = 10000)
sc.nor <- FindVariableFeatures(sc.nor, selection.method = "vst")
sc.nor <- ScaleData(sc.nor,features = rownames(sc.nor))


sc2 <- RunPCA(sc.nor,features = VariableFeatures(object = sc.nor), npcs = 100, verbose = FALSE)
DimPlot(sc2 ,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(sc2)
sc2 <- RunHarmony(sc2,"orig.ident", plot_convergence = TRUE)
sc2 <- FindNeighbors(sc2,dims = 1:50)
sc2 <- FindClusters(sc2,resolution = 0.5)
sc2 <- RunTSNE(sc2, dims = 1:50)
sc2 <- RunUMAP(sc2, dims = 1:50)
DimPlot(sc2, reduction = "umap", group.by = "ident", pt.size=0.5)
DimPlot(sc2, reduction = "umap", group.by = "ident", pt.size=0.5,split.by = 'orig.ident')


paramSweep_v2 <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}


#find best pK
sweep.res.list <- paramSweep_v2(sc2,PCs = 1:10,sct = F)
sweep.stats <- summarizeSweep(sweep.res.list,GT = F)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))


#ratio of doublet
annotations <- sc2$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc2@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))


doubletFinder_v2 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}

sc2_db <- doubletFinder_v2(sc2, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
sc2_db@meta.data[,"DF_hi.lo"] <- sc2_db@meta.data$DF.classifications_0.25_0.05_400
DimPlot(sc2_db,reduction = "umap",group.by = "DF_hi.lo")

sc3 <- subset(sc2_db, subset = DF_hi.lo == "Singlet")
VlnPlot(sc3, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,group.by = "orig.ident")
FeatureScatter(sc3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident")

sc3.nor <- NormalizeData(sc3, normalization.method = "LogNormalize", scale.factor = 10000)
sc3.nor <- FindVariableFeatures(sc3.nor, selection.method = "vst")
sc3.nor <- ScaleData(sc3.nor,features = rownames(sc.nor))

sc4 <- RunPCA(sc3.nor,features = VariableFeatures(object = sc3.nor), npcs = 100, verbose = FALSE)
DimPlot(sc4 ,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(sc4)
sc4 <- FindNeighbors(sc4,dims = 1:50)
sc4 <- FindClusters(sc4,resolution = 0.5)
sc4 <- RunTSNE(sc4, dims = 1:50)
sc4 <- RunUMAP(sc4, dims = 1:50)
DimPlot(sc4, reduction = "umap", group.by = "ident", pt.size=0.5)

CS_scRNA <- sc4
saveRDS(CS_scRNA,"CS_mergereps.rds")