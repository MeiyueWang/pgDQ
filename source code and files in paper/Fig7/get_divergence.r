library(rdist)
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
i <- "Development"
dfr <- fread(paste(i,"tpm.tsv",sep = "_"),stringsAsFactors = F,header = T,data.table = F)
rownames(dfr) <- dfr$gene
dfr <- dfr[,-1]
dfr.A <- dfr[orth$cs_A,];dfr.A[is.na(dfr.A)] <- 0;rownames(dfr.A) <- orth$groups
dfr.B <- dfr[orth$cs_B,];dfr.B[is.na(dfr.B)] <- 0;rownames(dfr.B) <- orth$groups
dfr.D <- dfr[orth$cs_D,];dfr.D[is.na(dfr.D)] <- 0;rownames(dfr.D) <- orth$groups
dfr.A <- dfr.A + 0.00001
dfr.B <- dfr.B + 0.00001
dfr.D <- dfr.D + 0.00001
dfr.A.1 <- dfr.A/(dfr.A+dfr.B+dfr.D)
dfr.B.1 <- dfr.B/(dfr.A+dfr.B+dfr.D)
dfr.D.1 <- dfr.D/(dfr.A+dfr.B+dfr.D)
df_allsubG <- rbind(dfr.A.1,dfr.B.1,dfr.D.1)
df_allsubG <- as.data.frame(df_allsubG)
list_allsubG <- as.list(df_allsubG)
list_allsubG_r <- lapply(list_allsubG,function(x){sum(as.numeric(rdist(t(matrix(x,nrow = nrow(orth),byrow = F)))))})
list_allsubG_r.df <- data.frame(samples=names(list_allsubG_r),divergence=as.numeric(list_allsubG_r))
write.table(list_allsubG_r.df,"divergence_bulk_tissue_development.txt",sep = "\t",row.names = F,col.names = T,quote = F)