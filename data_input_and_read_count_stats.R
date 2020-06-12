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

allowWGCNAThreads()

# Import kallisto files with txi 
# ==================================================================================
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
# ==================================================================================


# if already have a txi object, load it with the metadata (colData)
# load("~/Documents/bmc/Data/txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

# ==================================================================================
# read count stats chink starts here

expressed_genes<-txi.kallisto.tsv$counts
expressed_genes<-as.data.frame(expressed_genes)
expressed_genes$GeneID<-row.names(expressed_genes)
expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),]
expressed_genes<-expressed_genes[,c(49, 1:48)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:49)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Sample")
sub<-all_wheat_genes[,c(9, 2, 3, 4)]
rep_wise<-spread(sub, key = Rep, value=TPM)
rep_wise$Sum<-rep_wise$`1` + rep_wise$`2` + rep_wise$`3`
rep_wise$test1<-ifelse(rep_wise$`1`>0.5, 1,0)
rep_wise$test2<-ifelse(rep_wise$`2`>0.5, 1,0)
rep_wise$test3<-ifelse(rep_wise$`3`>0.5, 1,0)

# check correlation of reps
cor<-as.matrix(rep_wise[,c(3,4,5)])
cor<-rcorr(cor)
corrplot(cor$r, type="lower", order="original",p.mat = cor$P, 
         sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 2, 
         tl.srt = 0, tl.offset = 1, method="color", addCoef.col = "white")

rep_wise$sum2<-rep_wise$test1 + rep_wise$test2 + rep_wise$test3
expressed<-rep_wise[(rep_wise$sum2 >=2),] 
length(unique(expressed$GeneID))
unique_tx<-as.data.frame(unique(expressed$GeneID))
colnames(unique_tx)<-"GeneID"
unique_tx$Gene<-unique_tx$GeneID
unique_gene<- unique_tx %>% 
  separate(Gene, c("Gene", "Transcript"))
unique_gene$Tx<-unique_gene$GeneID
unique_gene<- unique_gene %>% 
  separate(Tx, c("dead", "Chromosome"), sep="S")
unique_gene$dead<-NULL
unique_gene <-unique_gene %>%
  separate(Chromosome, into=c("Chrom", "dead"), sep="0", extra="drop")
unique_gene$dead<-NULL
unique_gene$Chrom<-str_replace(unique_gene$Chrom, "U", "UU")
chroms<-unique_gene$Chrom
test<-strsplit(chroms, split="")
chroms<-do.call(rbind.data.frame, test)
colnames(chroms)<-c("Chromosome", "Genome")
unique_gene<-cbind(unique_gene, chroms)
write.csv(unique_gene, file="~/Documents/bmc/Data/all_expressed_genes.csv")
tab<-
  as.data.frame(table(unique_gene$Chromosome, unique_gene$Genome))
colnames(tab)<-c("Chromosome", "Genome", "Frequency")

ggplot(tab, aes(x=Chromosome, y=Frequency)) + 
  geom_col(aes(fill=Genome)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black")) +
  labs(fill = "Genome") +
  xlab("Chromosome") +
  ylab("Number of expressed genes / gene isoforms")+
  scale_fill_manual(values=c("grey80", "grey60", "grey40", "grey20"))


table<-merge(as.data.frame(table(expressed$Factor)),colData, by.x="Var1", by.y="Factor") 

ggplot(table, aes(x=as.factor(Timepoint), y=Freq)) +
  geom_col(aes(fill=Treatment), position="dodge") +
  theme_classic() +
  facet_grid(Genotype~.) +
  theme(legend.position = "right") + 
  scale_fill_manual(values=c("grey80","grey40"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Number of expressed genes / gene isoforms") + 
  xlab("Timepoint (dpi)")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black")) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(aes(yintercept=0)) 


# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype +  Rep + Treatment:Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Timepoint", "Genotype", "Rep"), returnData=TRUE)
pcaData$Sample<-paste(pcaData$Genotype, pcaData$Treatment)
# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2, group=Sample)) + geom_point(aes(colour=Sample), size= 7, alpha=0.9) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  scale_colour_manual(values=c("slategray3","dodgerblue", "rosybrown3", "orangered"), 
                      labels=c("Longbow + Tween20", expression(paste("Longbow + ", italic("Z. tritici"))), "Stigg + Tween20", expression(paste("Stigg + ", italic("Z. tritici"))))) +
  xlab("Principle component 1")+
  ylab("Principle component 2")


 # heatmap of genes by sample
dds<-estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("Genotype", "Treatment")])
pheatmap(assay(dds)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# bheatmap of sample distances
vst_counts<-assay(vsd)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Rep, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
