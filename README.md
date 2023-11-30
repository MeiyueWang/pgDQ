# pgDQ
 Pseudo-genome divergence quantification for polyploid wheat

## Pre-requisite
1) Seurat (v4.3.0)
2) rdist (v0.0.5)
3) Monocle (v2.22.0)

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

5) meristem marker genes as follows:
   | gene1 |
   | :---  |
   | gene2 | 
   | gene3 |
   | ...   |

## Run pgDQ

1) calculate subgenome divergence for all cells.
```js
Rscript 02divergence_cells.r /path/to/rds /path/to/orth
```
2) select initial state and calculate pseudotime for all cells.
Firstly, we choose cell clusters with highest expression values of meristem marker genes as initial clusters.
```js
Rscript 03initial_cluster.r /path/to/rds /path/to/marker_genes
```
