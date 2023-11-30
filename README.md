# pgDQ
 Pseudo-genome divergence quantification for polyploid wheat

## Pre-requisite
1) Seurat (v4.3.0)
2) 

## Input
1) run '01cluster_cells.r' to cluster cells and generate '.rds' data if you use 'filtered_feature_bc_matrix' from Cell Ranger output as input.
2) scRNA data in '.rds' format.
3) 1v1 orthologous genes across subgenomes as follows:
   | orthgroups | subgenome1 | subgenome2 | ... |
   | :--------- | :--------- | :--------- | :-- |
