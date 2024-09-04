library(dplyr)
library(ggplot2)
library(ggrepel)
dc <- read.table("ABD_distance_vs_pt_corr_rm_one_gene_tissue.txt")
dc <- dc[order(dc$V3,decreasing = F),]
dc$no <- c(1:nrow(dc))
dc$nor <- log(2/(1+dc$V3))
dc$nor_Z <- scale(dc$nor,scale = T,center = T)
wheat_func <- read.table("orth1v1_wheat_fucgene.txt",sep = "\t",stringsAsFactors = F)
rice_func <- read.table("orth1v1_rice_fucgene.txt",sep = "\t",stringsAsFactors = F)
og.sel.j <- left_join(dc,wheat_func,by=c("V1"="V1"))
og.sel.j <- left_join(og.sel.j,rice_func,by=c("V1"="V1"))
og.sel.j$func <- ifelse(is.na(og.sel.j$V2),NA,og.sel.j$V2)
og.sel.j$func <- ifelse(is.na(og.sel.j$V2.y),og.sel.j$func,og.sel.j$V2.y)
og.sel.j$func <- gsub(",.*","",og.sel.j$func)
orth <- read.table("Orthogroups_CS_ABD_1v1.txt",sep = "\t",stringsAsFactors = F,header = T)
orth$cs_A <- gsub("\\.[0-9]+","",orth$cs_A)
orth$cs_B <- gsub("\\.[0-9]+","",orth$cs_B)
orth$cs_D <- gsub("\\.[0-9]+","",orth$cs_D)
rownames(orth) <- orth$groups
p <- ggplot(og.sel.j)+geom_line(aes(no,nor_Z),color="grey")+
  geom_point(data=na.omit(og.sel.j[c(1:100),c("no","nor_Z","func")]),aes(no,nor_Z),color="black")+
  geom_text_repel(data=na.omit(og.sel.j[c(1:100),c("no","nor_Z","func")]),aes(no,nor_Z,label=func),max.overlaps = 35,color="black")+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Ranked gene list",y="",title = "")
ggsave("func_of_gene_by_diagnosis_tissue.pdf",p,height = 4.5,width = 4.5,units = "in",dpi = 300)