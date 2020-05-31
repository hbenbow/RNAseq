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
all_significant<- all_significant[- grep("LC", all_significant$ID),]


# lets check that all the DEGs are in the expressed gene set, and remove ones that are not
not_expressed_DEGs<-
  subset(all_significant, 
         !(all_significant$gene  %in% all_expressed_genes$GeneID))
all_significant<-all_significant[!(all_significant$gene==unique(not_expressed_DEGs$gene)),]
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
for(cultivar in c("Stigg", "Longbow")){
  # for(t in c(6, 24, 48, 96)){
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
# read in the output of the OmicsBox GO enrichment analysis and make graphs of the enriched categories
l48u <- read.delim("~/Documents/bmc/Data/DE_genes/l48u.txt")
l24u <- read.delim("~/Documents/bmc/Data/DE_genes/l24u.txt")

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
  theme(text = element_text(size=15, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 

# mp24_graph<-
ggplot(mp_24, aes(x = reorder(GO.Name, Difference), y=Percentage)) +
  geom_col(aes(fill=Set), position="dodge") + coord_flip() +
  scale_fill_manual(values=c("grey60", "grey40"),
                    labels=c("Reference Set", "SSP Set")) +
  theme_bw()+
  theme(text = element_text(size=15, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 

mp_48<-l48up[(l48up$GO.Category=="MOLECULAR_FUNCTION"),]
mp_48<-mp_48[order(-mp_48$Difference),]

# mp48_graph<-
ggplot(mp_48, aes(x = reorder(GO.Name, Difference), y=Percentage)) +
  geom_col(aes(fill=Set), position="dodge") + coord_flip() +
  scale_fill_manual(values=c("grey60", "grey40"),
                    labels=c("Reference Set", "SSP Set")) +
  theme_bw()+
  theme(text = element_text(size=15, colour="black"), 
        axis.text.x = element_text(colour="black"), 
        axis.text.y=element_text(colour="black"),
        legend.position = "none")+
  scale_x_discrete(labels = wrap_format(40))+
  scale_y_continuous()+
  xlab(("Gene Ontology Term"))+
  labs(fill = "Gene Set") 

# ==================================================================================
# clustering
# heatmap of genes by sample


# Import kallisto files with txi 

# Use these steps to import kallisto files
samples<-dir("~/Documents/bmc/Data/samples/") # where the directory 'samples'
# contains the kallisto output directories - 1 per sample.
files <- file.path(samples, "abundance.h5")
setwd("~/Documents/bmc/Data/samples/")
names(files) <- paste0(samples)
all(file.exists(files)) # need to navigate into samples directory for this to work!!
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.csv("~/Documents/bmc/Data/colData.csv")
setwd("~/Documents/bmc/Data/")
colData$Timepoint<-as.factor(colData$Timepoint)
colData$Rep<-as.factor(colData$Rep)

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype +  Rep + Treatment:Genotype)
dds<-estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("Genotype", "Treatment")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

