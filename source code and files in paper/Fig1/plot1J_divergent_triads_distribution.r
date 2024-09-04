library(ggplot2)
data <- read.table("cell_class_all_count.xls",sep = "\t",stringsAsFactors = F,header = T)
data$unb <- (16427-data$Balanced)/16427
p <- ggplot(data,aes(unb))+geom_density(color="#82AAE3",fill="#82AAE3",alpha = 0.5,size=1)+
  geom_vline(xintercept = 0.408413,linetype="dashed",size=1.5,color="#FEA1BF")+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(angle = 0,hjust = 0.5,vjust = 0.5))+
  labs(x="ratio of divergent triads",y="Density")
ggsave("distribution_of_divergent_triads_in_cells.pdf",p,width = 3.3,height = 3.3,dpi = 600)