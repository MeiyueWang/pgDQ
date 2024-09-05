library(ggalluvial)
library(dplyr)
library(openxlsx)
library(Seurat)
CS_scRNA <- readRDS("CS_mergereps.rds")
markers <- FindAllMarkers(CS_scRNA,only.pos = T,logfc.threshold=0.25)
write.table(markers,"markergenes_wheat.xls",sep = "\t",row.names = F,col.names = T,quote = F)

####
markers <- read.table("markergenes_wheat.xls",sep = "\t",stringsAsFactors = F,header = T)
markers <- markers[grep("TraesCS",markers$gene),]
markers_sel <- markers[markers$pct.1>0.25,]

no <- c()
for (i in unique(markers$cluster)) {
  no <- c(no,c(1:nrow(markers[markers$cluster==i,])))
}
markers$no <- no
#markers.sel <- markers[markers$no<=100,]
markers.sel <- markers_sel
wheat_pub <- read.xlsx("Asymmetric gene expression and cell-type-specific regulatory networks in the root of bread wheat revealed by single-cell multiomics analysis.xlsx",sheet = 2,colNames = F)
df_out <- c()
for (cus in c(0:19)) {
  cus.g <- markers.sel[markers.sel$cluster==cus,"gene"]
  for (pub in unique(sort(wheat_pub$X1))) {
    pub.g <- wheat_pub[wheat_pub$X1==pub,"X2"]
    a <- length(intersect(cus.g,pub.g))
    b <- length(cus.g) - a
    c <- length(pub.g) - a
    d <- 107891 - a - b - c
    fc <- (a/b)/(c/d)
    r <- c(cus,pub,length(cus.g),length(pub.g),fc)
    df_out <- rbind(df_out,r)
  }
}
df_out <- as.data.frame(df_out)
df_out$V5 <- as.numeric(df_out$V5)
mat <- matrix(df_out$V5,ncol = 22,byrow = T)
colnames(mat) <- paste("X",c(0:21),sep = "")
mat <- as.data.frame(mat)
mat$max <- apply(mat[,c(1:22)],1,function(x){names(x)[which(x==max(x))]})
mat$cluster <- c(0:19)

####AK58
cluster.AK58 <- as.character(c(0:21))
celltype.AK58 <- c("epidermis/cortex I","immatuer pericycle cells (IPC)","meristem I","proximal meristem","xylem pole pericycle (XPP)","provascular cells","root hair","meristem II","epidermis/root hair","phloem pole pericycle (PPP)","endoermis I (casparian strip)","metaxylem","epidermis/cortex II","protoxylem","immature sieve elements","root cap","protophloem","endodermis II (casparian strip)","stem cell niche (SCN)","root border cell","columella","companion cell")
df.AK58 <- data.frame(cluster.AK58,celltype.AK58,stringsAsFactors = F)
rownames(df.AK58) <- df.AK58$cluster.AK58
marker.AK58 <- read.xlsx("Asymmetric gene expression and cell-type-specific regulatory networks in the root of bread wheat revealed by single-cell multiomics analysis.xlsx",sheet = 2,colNames = F)
colnames(marker.AK58) <- c("cluster","mkgene")

####CS
cluster.CS <- as.character(c(0:19))
celltype.CS <- df.AK58[as.character(c(0,2,11,19,19,1,2,3,6,3,4,4,17,12,15,4,11,13,7,18)),"celltype.AK58"]
df.CS <- data.frame(cluster.CS,celltype.CS,stringsAsFactors = F)
marker.CS <- read.table("markergenes_wheat.xls",sep = "\t",stringsAsFactors = F,header = T)
marker.CS <- marker.CS[grep("TraesCS",marker.CS$gene),]
marker.CS <- marker.CS[marker.CS$pct.1>0.25,]
marker.CS <- marker.CS[,c("cluster","gene")]
colnames(marker.CS) <- c("cluster","mkgene")


out.AK58 <- c()
for (i in unique(sort(df.CS$celltype.CS))) {
  mkgenes.i <- marker.CS[marker.CS$cluster%in%df.CS[df.CS$celltype.CS==i,"cluster.CS"],"mkgene"]
  for (j in unique(sort(df.AK58$celltype.AK58))) {
    mkgenes.j <- marker.AK58[marker.AK58$cluster%in%df.AK58[df.AK58$celltype.AK58==j,"cluster.AK58"],"mkgene"]
    Freq.ij <- length(intersect(mkgenes.i,mkgenes.j))
    r <- c(i,j,Freq.ij)
    out.AK58 <- rbind(out.AK58,r)
  }
}
out.AK58 <- as.data.frame(out.AK58,stringsAsFactors = F)
out.AK58$V3 <- as.numeric(out.AK58$V3)
colnames(out.AK58) <- c("CS","AK58","Freq")
p <- ggplot(out.AK58,aes(y=Freq,axis1=CS,axis2=AK58))+
  geom_alluvium(aes(fill=CS))+
  geom_stratum(width = 1/12) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("CS", "AK58"), expand = c(.05, .05)) +
  theme_bw()


ggsave("celltype_annotation.pdf",p,width = 7.38,height = 6.8,dpi = 300)
