library(Seurat)
library(monocle)

args <- commandArgs(T)
root_state <- as.numeric(args[1])

load("pgDQ_monocle.RData")
cds <- orderCells(cds,root_state = root_state)
cds.data <- pData(cds)
cds.data$cells <- rownames(cds.data)
write.table(cds.data,"pgDQ_monocle_pseudotime.xls",sep = "\t",quote = F,row.names = F,col.names = T)
