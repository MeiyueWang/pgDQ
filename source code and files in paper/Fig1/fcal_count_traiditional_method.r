library(data.table)
fcount <- function(x){
  x <- as.character(x["class"])
  x <- unlist(strsplit(x,","))
  n <- data.frame(table(x),stringsAsFactors = F)
  rownames(n) <- n$x
  n <- n[c("1A dominant","2A suppressed","3B dominant","4B suppressed","5D dominant","6D suppressed","7Balanced"),]
  r <- paste0(n$Freq,collapse = ",")
  return(r)
}

data <- fread("inr_0_unbalanced_cell_class_all.xls",sep = "\t",header = T,data.table = F)
data$count <- apply(data,1,fcount)
write.table(data,"cell_class_all_count.xls",row.names = F,col.names = T,quote = F,sep = "\t")
