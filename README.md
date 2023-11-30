# pgDQ
 Pseudo-genome divergence quantification for polyploid wheat

## Pre-requisite
1) Seurat (v4.3.0)
2) 

## Input
1) run '01cluster_cells.r' to cluster cells and generate '.rds' data if you use 'filtered_feature_bc_matrix' from Cell Ranger output as input.
   ```js
   Rscript 01cluster_cells.r /path/to/filtered_feature_bc_matrix
   ```
3) scRNA data in '.rds' format.
4) 1v1 orthologous genes across subgenomes as follows:
   
   | orthgroups | subgenomeA | subgenomeB | subgenomeC | ... |
   | :--------- | :--------- | :--------- | :--------- | :-- |
   | orth1      | geneA1     | geneB1     | geneC1     | ... |
   | orth2      | geneA2     | geneB2     | geneC2     | ... |

## Run pgDQ

