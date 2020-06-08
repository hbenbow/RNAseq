all_significant <- 
  read.csv("~/Documents/bmc/Data/DE_genes/all_significant.csv", row.names=1)
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
all_significant$Factor<-paste(all_significant$Cultivar, all_significant$Timepoint, all_significant$Regulation)

list<-list()
for(factor in unique(all_significant$Factor)){
  data<-all_significant[all_significant$Factor==factor,]
  data$abs<-abs(data$log2FoldChange)
  top5<-top_n(data, 5, data$abs)
  list[[(length(list)+1)]]<-top5
}
all<-as.data.frame(do.call("rbind", list))

