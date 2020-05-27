all_significant <- 
  read.csv("~/Documents/bmc/Data/all_significant.csv", row.names=1)
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
all_significant<- all_significant[- grep("LC", all_significant$ID),]
