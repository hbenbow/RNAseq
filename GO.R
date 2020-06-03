setwd("~/Documents/bmc/Data/DE_genes/B2G/")
list<-list()
for(i in dir(pattern="txt")){
  data<-read.delim(i, row.names=NULL)
  name<-paste(i)
  name<-substr(i, 23, (nchar(i)-4))
  things<-strsplit(name, "_")[[1]]
  data$Timepoint<-things[1]
  data$regulation<-things[2]
  data$Cultivar<-things[3]
  data$Cultivar<-sub("L", "Longbow", data$Cultivar)
  data$Cultivar<-sub("S", "Stigg", data$Cultivar)
  assign(name, data)
  rm(i)
  rm(data)
  rm(name)
  rm(things)
  list[[(length(list)+1)]]<-data
}

all<-do.call(rbind.data.frame(list))
