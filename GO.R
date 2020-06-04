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
options(scipen = 999)

allowWGCNAThreads()



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
  rm(name)
  rm(things)
  list[[length(list)+1]]<-data
  rm(data)
}

all<-do.call(rbind.data.frame, list)
write.csv(all, file="~/Documents/bmc/Data/DE_genes/B2G/all_bp.csv")
# at this point, open the file manually, and assign high level processes to each BP

all_bp <- read.csv("~/Documents/bmc/Data/DE_genes/B2G/all_bp.csv")
all_bp$regulation<-sub("up", "Up", all_bp$regulation)
all_bp$regulation<-sub("down", "Down", all_bp$regulation)

bp_for_plot<-aggregate(all_bp$GO, by=list(all_bp$H_L_process, all_bp$Cultivar, all_bp$regulation), FUN="sum")
colnames(bp_for_plot)<-c("H_L_process", "Cultivar", "regulation", "No.Seqs")

levels<-aggregate(bp_for_plot$No.Seqs, by=list(bp_for_plot$H_L_process), FUN="sum")
levels<-levels[order(levels$Group.1),]
bp_for_plot$H_L_process <- factor(bp_for_plot$H_L_process, levels = c("Other", "Growth/development", "Biosynthesis",
                                                                      "Redox", "Catabolism",
                                                                      "Transport", "Hormone",
                                                                      "PTM", "Metabolism", "Stress reponse",  
                                                                      "Photosynthesis", "Transcription"))


ggplot(bp_for_plot, aes(x=H_L_process, y=No.Seqs, fill=Cultivar)) +
  geom_col(aes(fill=factor(Cultivar, levels=c("Longbow", "Stigg"))), position=position_dodge(preserve="single")) +
  facet_wrap(regulation~.) +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=20, colour='black'), axis.text.x = element_text(colour="black")) +
  scale_fill_manual(values=c("grey30", "grey60"), limits=c("Longbow", "Stigg")) +
  ylab("Number of differentially expressed genes") +
  xlab("High level biological process")


other<-all_bp[(all_bp$H_L_process=="Other"),]
other_Stigg<-other[(other$Cultivar=="Stigg"),]
other_Stigg<-as.character(other_Stigg$GO_process)
other_Longbow<-other[(other$Cultivar=="Longbow"),]
other_Longbow<-as.character(other_Longbow$GO_process)

stigg_BP_other<-subset(other_Stigg, !(other_Stigg %in% other_Longbow))
stigg_BP_other<-subset(other, other$GO_process %in% stigg_BP_other)

ggplot(stigg_BP_other, aes(x=GO_process, y=GO, fill=regulation)) +
  geom_col(aes(fill=factor(regulation, levels=c("Down", "Up"))), position=position_dodge(preserve="single")) +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=15, colour='black'), axis.text.x = element_text(colour="black")) +
  scale_fill_manual(values=c("grey30", "grey60"), limits=c("Down", "Up")) +
  scale_x_discrete(labels=wrap_format(40))+
  ylab("Number of differentially expressed genes") +
  xlab("biological process")

write.csv(stigg_BP_other, file="~/Documents/bmc/Data/DE_genes/B2G/Stigg_bp_other.csv")

stigg_specific_up<-read.delim("~/Documents/bmc/Data/DE_genes/Stigg_sepcific/up_bp.txt")
stigg_specific_down<-read.delim("~/Documents/bmc/Data/DE_genes/Stigg_sepcific/Down_bp.txt")
Longbow_specific_up<-read.delim("~/Documents/bmc/Data/DE_genes/Longbow_specific/Up_bp.txt")
Longbow_specific_down<-read.delim("~/Documents/bmc/Data/DE_genes/Longbow_specific/Down_bp.txt")
stigg_specific_down$Go_term<-row.names(stigg_specific_down)
stigg_specific_up$Go_term<-row.names(stigg_specific_up)
Longbow_specific_down$Go_term<-row.names(Longbow_specific_down)
Longbow_specific_up$Go_term<-row.names(Longbow_specific_up)


ggplot(stigg_specific_down, aes(x=reorder(Go_term, GO), y=GO)) +
  geom_col() +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=17, colour='black'), axis.text.x = element_text(colour="black"))+
  scale_x_discrete(labels=wrap_format(50))+
  ylab("Number of differentially expressed genes") +
  xlab("Biological process")

ggplot(stigg_specific_up, aes(x=reorder(Go_term, GO), y=GO)) +
  geom_col() +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=17, colour='black'), axis.text.x = element_text(colour="black"))+
  scale_x_discrete(labels=wrap_format(50))+
  ylab("Number of differentially expressed genes") +
  xlab("Biological process")
