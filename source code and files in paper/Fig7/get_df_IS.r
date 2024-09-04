dc <- read.table("ABD_distance_vs_pt_corr_rm_one_gene_tissue.txt")
dc <- dc[order(dc$V3,decreasing = F),]
dc$nor <- log(2/(1+dc$V3))
dc$nor_Z <- scale(dc$nor,scale = T,center = T)
dc$no <- c(1:nrow(dc))
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
dc.A <- dc
dc.A$NAME <- orth[dc$V1,"cs_A"]
dc.B <- dc
dc.B$NAME <- orth[dc$V1,"cs_B"]
dc.D <- dc
dc.D$NAME <- orth[dc$V1,"cs_D"]
dc_out <- rbind(dc.A,dc.B,dc.D)
dc_out$DESCRIPTION <- "na"
dc_out$zero <- 0
dc_out <- dc_out[order(dc_out$no),]
dc_out <- dc_out[,c("NAME","DESCRIPTION","nor_Z","zero")]
colnames(dc_out) <- c("NAME","DESCRIPTION","IS","zero")
write.table(dc_out,"df_IS_tissue.txt",sep = "\t",col.names = T,row.names = F,quote = F)