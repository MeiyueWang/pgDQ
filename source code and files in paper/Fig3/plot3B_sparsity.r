library(ggplot2)
dfpg <- c()
dfdtb <- c()
for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) {
  pg <- read.table(paste("sparsity_simu150/simu_inr",i,"pgDQ_corr_sparsity.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  pg$percent <- i
  dfpg <- rbind(dfpg,pg)
  dtb <- read.table(paste("sparsity_simu150/simu_inr",i,"dtb_corr_sparsity.xls",sep = "_"),sep = "\t",stringsAsFactors = F)
  dtb$percent <- i
  dfdtb <- rbind(dfdtb,dtb)
}
dfpg$method <- "pgDQ"
dfdtb$method <- "dtb"
df <- rbind(dfpg,dfdtb)
df.ag <- aggregate(.~percent+method,df,mean)
df.ag$method <- factor(df.ag$method,levels = c("pgDQ","dtb"))
p <- ggplot(df.ag,aes(percent,V1,color=method))+geom_point(size=5)+
  geom_line(size=1)+scale_color_manual(values = c("#A51E26","#162A77"))+
  scale_x_continuous(breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                     labels = as.character(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)))+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(color="",x="",y="")+ylim(0.2,1)
ggsave("plot3B_evaluation_simulated_sparsity.pdf",p,height = 3.6,width = 5.1,units = "in",dpi = 300)
