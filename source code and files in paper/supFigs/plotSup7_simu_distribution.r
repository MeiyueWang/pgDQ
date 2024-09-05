library(data.table)
library(ggplot2)
setwd("D:/05.scRNA-seq_jpgDR/TableCollection/sparsity_simu150")
dfpg <- c()
dfdtb <- c()
for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) {
  pg <- read.table(paste("simu_inr",i,"pgDQ_corr_sparsity.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  pg$percent <- i
  dfpg <- rbind(dfpg,pg)
  dtb <- read.table(paste("simu_inr",i,"dtb_corr_sparsity.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  dtb$percent <- i
  dfdtb <- rbind(dfdtb,dtb)
}
dfpg$method <- "pgDQ"
dfdtb$method <- "dtb"
df <- rbind(dfpg,dfdtb)
df$method <- factor(df$method,levels = c("pgDQ","dtb"))
p1 <- ggplot(df,aes(V1,color=method))+geom_freqpoly()+
  facet_grid(rows = vars(percent))+
  scale_color_manual(values = c("#A51E26","#162A77"))+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(colour = "black",size = 12),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",x="",y="")+guides(color="none")
ggsave("distribution_of_sparsity_simu150.pdf",p1,width = 4,height = 6,dpi = 300)


setwd("noise_simu150")
dfpg <- c()
dfdtb <- c()
for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) {
  pg <- read.table(paste("simu_inr",i,"pgDQ_corr_noise.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  pg$percent <- i
  dfpg <- rbind(dfpg,pg)
  dtb <- read.table(paste("simu_inr",i,"dtb_corr_noise.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  dtb$percent <- i
  dfdtb <- rbind(dfdtb,dtb)
}
dfpg$method <- "pgDQ"
dfdtb$method <- "dtb"
df <- rbind(dfpg,dfdtb)
df$method <- factor(df$method,levels = c("pgDQ","dtb"))
p2 <- ggplot(df,aes(V1,color=method))+geom_freqpoly()+
  facet_grid(rows = vars(percent))+
  scale_color_manual(values = c("#A51E26","#162A77"))+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(colour = "black",size = 12),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(color="",x="",y="")
ggsave("distribution_of_noise_simu150.pdf",p2,width = 4,height = 6,dpi = 300)