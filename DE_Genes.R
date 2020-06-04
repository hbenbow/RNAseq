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
library(Hmisc)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(scales)
options(scipen = 999)

allowWGCNAThreads()

# all DE genes from larger DE analysis set, extract subgroups of DE genes by timepoint and cultivar

all_expressed_genes <- read.csv("~/Documents/bmc/Data/all_expressed_genes.csv", row.names=1)
all_significant <- 
  read.csv("~/Documents/bmc/Data/DE_genes/all_significant.csv", row.names=1)
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
# all_significant<- all_significant[- grep("LC", all_significant$ID),]


# lets check that all the DEGs are in the expressed gene set, and remove ones that are not
not_expressed_DEGs<-
  subset(all_significant, 
         !(all_significant$gene  %in% all_expressed_genes$GeneID))
# all_significant<-all_significant[!(all_significant$gene==unique(not_expressed_DEGs$gene)),]
# ==================================================================
# lets look at genome distribution of DE genes
all_significant$Gene<-all_significant$gene
all_significant<- all_significant %>% 
  separate(Gene, c("dead", "Chromosome"), sep="S")
all_significant$dead<-NULL
all_significant <-all_significant %>%
  separate(Chromosome, into=c("Chrom", "dead"), sep="0", extra="drop")
all_significant$dead<-NULL
all_significant$Chrom<-str_replace(all_significant$Chrom, "U", "UU")
chroms<-all_significant$Chrom
test<-strsplit(chroms, split="")
chroms<-do.call(rbind.data.frame, test)
colnames(chroms)<-c("Chromosome", "Genome")
all_significant<-cbind(all_significant, chroms)

tab<-
  as.data.frame(table(all_significant$Chromosome, all_significant$Genome))
colnames(tab)<-c("Chromosome", "Genome", "Frequency")

write.csv(all_significant, file="~/Documents/bmc/Data/DE_genes/all_significant.csv")

# figure 3 sort this out
ggplot(genome_percentage, aes(x=Cultivar, y=Percentage, group=Genome)) + 
  theme_bw() +
  geom_bar(stat="identity", aes(fill=Genome), position="dodge") + 
  scale_fill_manual(values=c("grey80", "grey60", "grey40", "black")) + ylab("Percentage of differentially expressed genes")  +
  facet_grid(Regulation~Timepoint) +
  theme(text = element_text(size=15, colour="black")) +
  geom_point(aes(y=Expected,group=Genome),colour="blue", size=3, position = position_dodge(width=0.9)) +
  geom_hline(aes(yintercept=32))

# ==================================================================
# export plain text files of lists of genes that are DE in each timepoint/cultivar/regulation combination
# these files can be used in OmicsBox for GO enrichment analysis
# for(cultivar in c("Stigg", "Longbow")){
for(t in c(6, 24, 48, 96)){
  for(r in c("Up", "Down")){
    data<-all_significant[(all_significant$Cultivar==cultivar),]
    data<-data[(data$Timepoint==t),]
    data<-data[(data$Regulation==r),]
    assign(paste(cultivar, t, r), data)
    genes<-data$ID
    write.table(genes, file=paste("~/Documents/bmc/Data/DE_genes/", cultivar, r, ".txt", sep=""), row.names = F, quote=F, col.names = F)
  }
}
# }

# ==================================================================
# subset genes that are DE in Stigg and not in Longbow
Stigg<-all_significant[(all_significant$Cultivar=="Stigg"),]
Longbow<-all_significant[(all_significant$Cultivar=="Longbow"),]
Stigg<-Stigg$ID
Longbow<-Longbow$ID

Unique_Stigg<-subset(Stigg, !(Stigg %in% Longbow))
Unique_Stigg_DEG<-subset(all_significant, all_significant$ID %in% Unique_Stigg)

table(Unique_Stigg_DEG$Chromosome, Unique_Stigg_DEG$Genome.3)

# for(t in c()){
for(r in c("Up", "Down")){
  data<-Unique_Stigg_DEG[(Unique_Stigg_DEG$Regulation==r),]
  # data<-data[(data$Regulation==r),]
  assign(paste(t, r), data)
  genes<-data$ID
  write.table(genes, file=paste("~/Documents/bmc/Data/DE_genes/Stigg_sepcific/",  r, ".txt", sep=""), row.names = F, quote=F, col.names = F)
}
# }

# ==================================================================
# subset genes that are DE in Stigg and not in Longbow
Stigg<-all_significant[(all_significant$Cultivar=="Stigg"),]
Longbow<-all_significant[(all_significant$Cultivar=="Longbow"),]
Stigg<-Stigg$ID
Longbow<-Longbow$ID

Unique_Longbow<-subset(Longbow, !(Longbow %in% Stigg))
Unique_Longbow_DEG<-subset(all_significant, all_significant$ID %in% Unique_Longbow)

table(Unique_Longbow_DEG$Chromosome, Unique_Longbow_DEG$Genome.3)

# for(t in c()){
for(r in c("Up", "Down")){
  data<-Unique_Longbow_DEG[(Unique_Longbow_DEG$Regulation==r),]
  # data<-data[(data$Regulation==r),]
  assign(paste(t, r), data)
  genes<-data$ID
  write.table(genes, file=paste("~/Documents/bmc/Data/DE_genes/Longbow_specific/",  r, ".txt", sep=""), row.names = F, quote=F, col.names = F)
}
# }

L_up<-read.delim("~/Documents/bmc/Data/DE_genes/Longbow_specific/Up_bp.txt")
L_up$Cultivar<-"Longbow"
L_up$Regulation<-"Up"
L_up$GO_Category<-row.names(L_up)
L_down<-read.delim("~/Documents/bmc/Data/DE_genes/Longbow_specific/Down_bp.txt")
L_down$Cultivar<-"Longbow"
L_down$Regulation<-"Down"
L_down$GO_Category<-row.names(L_down)
S_up<-read.delim("~/Documents/bmc/Data/DE_genes/Stigg_sepcific/up_bp.txt")
S_up$Cultivar<-"Stigg"
S_up$Regulation<-"Up"
S_up$GO_Category<-row.names(S_up)
S_down<-read.delim("~/Documents/bmc/Data/DE_genes/Stigg_sepcific/down_bp.txt")
S_down$Cultivar<-"Stigg"
S_down$Regulation<-"Down"
S_down$GO_Category<-row.names(S_down)


all_bp<-rbind(L_up, L_down, S_up, S_down)
all_up<-all_bp[(all_bp$Regulation=="Up"),]
all_down<-all_bp[(all_bp$Regulation=="Down"),]

write.csv(all_bp, file="~/Documents/bmc/Data/DE_genes/all_bp.csv")
#   plot<-
ggplot(S_up[(S_up$GO>=2),], aes(x=as.factor(GO_Category), y=GO)) + 
  geom_col(aes(fill=Cultivar), position = position_dodge(preserve = "single")) + 
  # facet_wrap(~Regulation, ncol=1) + 
  theme_classic() +
  scale_fill_manual(values=c("grey30", "grey50")) +
  xlab("Gene Ontology Term") +
  theme(text = element_text(size=15, colour='black'), axis.text.x = element_text(colour="black")) +
  ylab("Number of differentially expressed genes") +
  scale_x_discrete(labels = wrap_format(40))+
  coord_flip()
# ggsave(plot, file=paste(i, ".pdf"))
ggplot(L_up[(L_up$GO>=3),], aes(x=as.factor(GO_Category), y=GO)) + 
  geom_col(aes(fill=Cultivar), position = position_dodge(preserve = "single")) + 
  # facet_wrap(~Regulation, ncol=1) + 
  theme_classic() +
  scale_fill_manual(values=c("grey30", "grey50")) +
  xlab("Gene Ontology Term") +
  theme(text = element_text(size=15, colour='black'), axis.text.x = element_text(colour="black")) +
  ylab("Number of differentially expressed genes") +
  scale_x_discrete(labels = wrap_format(40))+
  coord_flip()

# ==================================================================
# read in the output of the OmicsBox GO enrichment analysis and make graphs of the enriched categories
l48u <- read.delim("~/Documents/bmc/Data/DE_genes/All_Subsets/l48u.txt")
l24u <- read.delim("~/Documents/bmc/Data/DE_genes/All_Subsets/l24u.txt")

l24u$Test<-l24u$Nr.Test/l24u$Non.Annot.Test*100
l24u$Reference<-l24u$Nr.Reference/l24u$Non.Annot.Reference*100
l24u$Difference<-l24u$Test-l24u$Reference
l24test<-l24u[,c(1:11,13)] 
l24test$Set<-"Test"
l24ref<-l24u[,c(1:10, 12,13)] 
l24ref$Set<-"Reference"
colnames(l24ref)<-colnames(l24test)
l24up<-rbind(l24ref, l24test)
colnames(l24up)[11]<-"Percentage"


l48u$Test<-l48u$Nr.Test/l48u$Non.Annot.Test*100
l48u$Reference<-l48u$Nr.Reference/l48u$Non.Annot.Reference*100
l48u$Difference<-l48u$Test-l48u$Reference
l48test<-l48u[,c(1:11,13)] 
l48test$Set<-"Test"
l48ref<-l48u[,c(1:10, 12,13)] 
l48ref$Set<-"Reference"
colnames(l48ref)<-colnames(l48test)
l48up<-rbind(l48ref, l48test)
colnames(l48up)[11]<-"Percentage"

bp_24<-l24up[(l24up$GO.Category=="BIOLOGICAL_PROCESS"),]
bp_24<-bp_24[order(-bp_24$Difference),]

mp_24<-l24up[(l24up$GO.Category=="MOLECULAR_FUNCTION"),]
mp_24<-mp_24[order(-mp_24$Difference),]

# bp24_graph<-
ggplot(bp_24[1:40,], aes(x = reorder(GO.Name, Difference), y=Percentage)) +
  geom_col(aes(fill=Set), position="dodge") + coord_flip() +
  scale_fill_manual(values=c("grey60", "grey40"),
                    labels=c("Reference Set", "SSP Set")) +
  theme_bw()+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 
ggsave(file="~/Documents/bmc/Graphs/bp24.pdf")

# mp24_graph<-
ggplot(mp_24, aes(x = reorder(GO.Name, Difference), y=Percentage)) +
  geom_col(aes(fill=Set), position="dodge") + coord_flip() +
  scale_fill_manual(values=c("grey60", "grey40"),
                    labels=c("Reference Set", "SSP Set")) +
  theme_bw()+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 
ggsave(file="~/Documents/bmc/Graphs/mf24.pdf")

mp_48<-l48up[(l48up$GO.Category=="MOLECULAR_FUNCTION"),]
mp_48<-mp_48[order(-mp_48$Difference),]

# mp48_graph<-
ggplot(mp_48, aes(x = reorder(GO.Name, Difference), y=Percentage)) +
  geom_col(aes(fill=Set), position="dodge") + coord_flip() +
  scale_fill_manual(values=c("grey60", "grey40"),
                    labels=c("Reference Set", "Test Set")) +
  theme_bw()+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 
ggsave(file="~/Documents/bmc/Graphs/legend.pdf")

# ==================================================================
# WRKY TF data

wrkys<-c("TraesCS3A02G347500.1", 
         "TraesCS3B02G379200.1",
         "TraesCS3D02G341100.1",
         "TraesCS6A02G146900.1",
         "TraesCS6B02G175100.2",
         "TraesCS6D02G130600.1")
wrkys_exp<-subset(all_significant, all_significant$ID %in% wrkys)

wrkys_exp$Gene_group<-as.factor(c(1, 1, 1, 2, 2, 3))

ggplot(wrkys_exp, aes(x=ID, y=log2FoldChange)) +
  geom_bar(stat="identity", aes(fill=Gene_group)) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), 
                width=.5) +
  theme_classic() +
  scale_fill_manual(values=c("grey30", "grey50", "grey70"),
                    labels=c("Group 1", "Group 2", "Group 3")) +
  theme(text = element_text(size=15, colour="black"), 
        axis.text.x = element_text(colour="black", angle = 50, hjust=1), 
        axis.text.y=element_text(colour="black")) + 
  labs(fill="Homologous group") +
  xlab("Gene ID") + ylab("Log(2) Fold Change")
  
