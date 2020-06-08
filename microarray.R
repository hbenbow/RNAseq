library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(dplyr)
library(plyr)
library(Hmisc)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(scales)
setwd("~/Documents/bmc/Data/microarray/")
Microarray <- read.csv("~/Documents/bmc/Data/microarray/Microarray dataset 11.06.19.csv")
blast <- read.delim("~/Documents/bmc/Data/microarray/probe_to_cdsv1.1.tab", header=F)
all_significant <- 
  read.csv("~/Documents/bmc/Data/DE_genes/all_significant.csv", row.names=1)
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
all_significant$Factor<-paste(all_significant$Cultivar, all_significant$Timepoint)
colnames(blast)<-c(
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",	
  "gapopen",
  "qstart",
  "qend",
  "sstart",	 
  "send",
  "evalue",
  "bitscore")	

attach(blast)
blast.sorted<-blast[order(qseqid, -bitscore, evalue, -pident),]
blast.sorted$hsp<-paste(blast.sorted$qseqid, blast.sorted$sseqid)
blast.sorted<-blast.sorted[!(duplicated(blast.sorted$hsp)),]
detach(blast)

top6<-list()
for(query in unique(blast.sorted$qseqid)){
  data<-blast.sorted[(blast.sorted$qseqid==query),]
  top<-top_n(data, 6, row.names(data))
  top6[[(length(top6) + 1)]]=top
}

top_6<-as.data.frame(do.call("rbind", top6))

top_6$qseqid<-gsub(":", "", top_6$qseqid)

list<-list()
for(i in unique(all_significant$ID)){
  data<-all_significant[all_significant$ID==i,]
  data<-data[,c(9,2,20)]
  colnames(data)<-c("ID", "LFC", "Condition")
  data$Data<-"RNAseq"
  gene<-as.character(data$ID)
  probes<-subset(top_6, top_6$sseqid %in% gene)
  if(length(probes)>1){
    for(prob in unique(probes$qseqid)){
      probe<-as.character(prob)
      mic_data<-Microarray[Microarray$ProbeID==probe,]
      if(nrow(mic_data)>0){
        mic_data<-mic_data[,c(1, 2, 10)]
        colnames(mic_data)<-c("ID", "LFC", "Condition")
        mic_data$Data<-"Microarray"
        d2<-rbind(data, mic_data)
        list[[(length(list) + 1)]]=d2
      }
    }
  }
}

for(i in 1:length(list)){
  list[[i]]$Pair<-paste(i)
}

validation<-as.data.frame(do.call("rbind", list))



write.csv(top_6, file="~/Documents/bmc/Data/microarray/top_6.csv")
genes<-top_3$sseqid
write.table(genes, file=paste(sets, "_genes.txt", sep=""), col.names=F, row.names = F, quote=F)

