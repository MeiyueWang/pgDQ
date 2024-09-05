sample <- c("rep1","rep2")
uniq <- c("314513269","233731940")
multi <- c("16326595","14336926")
df <- data.frame(sample=sample,uniq=uniq,multi=multi,stringsAsFactors = F)
df_melt <- reshape2::melt(df,id.vars=c("sample"),measure.vars=c("uniq","multi"))
df_melt$value <- as.numeric(df_melt$value)
p <- ggplot(df_melt,aes(sample,value,fill=variable))+geom_bar(stat = "identity")+
  theme_classic()+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(fill="",title = "",x="Sample",y="# mapped reads")
ggsave("cellranger_mapping_stats.pdf",p,width = 3.9,height = 3.6,dpi = 300)


314513269/(314513269+16326595)
16326595/(314513269+16326595)


233731940/(233731940+14336926)
14336926/(233731940+14336926)