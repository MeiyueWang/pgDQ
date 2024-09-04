library(dplyr)
library(ggplot2)
library(ggrepel)
dc <- read.table("divGene_out.20240412.xls",sep = "\t",stringsAsFactors = F)
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
og.sel.j[og.sel.j$V1=="OG0014846","func"] <- "TaTIP1"
p <- ggplot(og.sel.j,aes(no,nor_Z,label=func))+geom_point(color="red")+geom_line()+
  geom_text_repel(max.overlaps = 35,color="black")+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="Ranked gene list",y="log10(1-corr)",title = "")
ggsave("funcGenes_of_divGene_out.pdf",p,height = 4.5,width = 4.5,units = "in",dpi = 300)