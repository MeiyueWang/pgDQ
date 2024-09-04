library(ggplot2)
df <- c()
for (i in c(1:1000)) {
  fdr <- 1-0.95^i
  r <- c(i,fdr)
  df <- rbind(df,r)
}
df <- as.data.frame(df,stringsAsFactors = F)
df$class <- ifelse(df$V1%in%c(10,50,100,1000),"sel","nonsel")
p <- ggplot()+geom_point(data=subset(df,class=="sel"),aes(V1,V2))+
  geom_line(data=df,aes(V1,V2),color="red",size=1)+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(colour = "black",size = 14),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1,colour = "black"),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))+
  labs(x="# triad genes",y="Family-wise error rate")
ggsave("error_rate_distribution.pdf",p,height = 5,width = 5)