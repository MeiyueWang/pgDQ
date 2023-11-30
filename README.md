# pgDQ
 Pseudo-genome divergence quantification for polyploid wheat

## Pre-requisite
1) Seurat (v4.3.0)
2) rdist (v0.0.5)
3) Monocle (v2.22.0)
4) foreach (v1.5.2)
5) doParallel (v1.0.17)

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
   Firstly, we choose cell cluster with highest expression values of meristem marker genes as initial cluster.
   ```js
   Rscript 03initial_cluster.r /path/to/rds /path/to/marker_genes
   Rscript 04plot_trajectory_monocle.r /path/to/rds
   ```
   then, get the initial cluster from output 'pgDQ_initial_cluster.txt', and find the root state from trajectory plot:
   ![pic_cluster](https://github.com/MeiyueWang/pgDQ/blob/main/trajectory_colored_by_clusters.png)
   ![pic_state](https://github.com/MeiyueWang/pgDQ/blob/main/trajectory_colored_by_states.png)

   Lastly, order cells with root state and calculate pseudotime for all cells.
   ```js
   Rscript 05pseudotime_of_cells.r root_state
   ```
3) divGene algorithm to calculate the indispensable scores for all orthologous gene groups.
   ```js
   Rscript 06divGene.r /path/to/rds /path/to/orth n_cores
   ```

## Output

1) output of '01cluster_cells.r' is 'pgDQ_scRNA.rds'.
   
2) output of '02divergence_of_cells.r' is 'pgDQ_subgenome_divergence_of_all_cells.xls':
   | V1    | V2   |
   | :---- | :--- |
   | cell1 | div1 | 
   | cell2 | div2 |
   | ...   | ...  |
   
3) output of '03initial_cluster.r' is 'pgDQ_initial_cluster.txt'.

4) output of '04plot_trajectory_monocle.r' is 'pgDQ_monocle.RData','pgDQ_monocle_trajectory_colored_by_cell_clusters.png' and 'pgDQ_monocle_trajectory_colored_by_state.png'.

5) output of '05pseudotime_of_cells.r' is 'pgDQ_monocle_pseudotime.xls'.

6) output of '06divGene.r' is 'pgDQ_indispensable_scores_of_all_1v1orthologous_genes.xls':
   | orthgroups    | subgenomeA | subgenomeB | subgenomeC | ... | pvalue | corr | IS |
   | :------------ | :--------- | :--------- | :--------- | :-- | :----- | :--- | :- |
   | orth1         | geneA1     | geneB1     | geneC1     | ... | p1     | corr1| IS1|   
   | orth2         | geneA2     | geneB2     | geneC2     | ... | p2     | corr2| IS2|
   | ...           | ...        | ...        | ...        | ... | ...    | ...  | ...|

  
   
