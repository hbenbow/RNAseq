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
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
options(scipen = 999)

allowWGCNAThreads()

# all DE genes from larger DE analysis set, extract subgroups of DE genes by timepoint and cultivar

all_expressed_genes <- read.csv("~/Documents/bmc/Data/all_expressed_genes.csv", row.names=1)
all_significant <- 
  read.csv(file="~/Documents/bmc/Data/DE_genes/all_significant.csv", row.names=1)
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]

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
rm(files)
rm(samples)

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~ Treatment + 
                                  Timepoint + Genotype + Treatment:Genotype)
vst<-varianceStabilizingTransformation(dds)
vds<-assay(vst)

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))
rm(order)
counts1<-txi.kallisto.tsv$abundance
counts1<-scale(counts1)
counts1<-vds
counts<-counts1
colData$TxGxT<-paste(colData$Rep, colData$Genotype, colData$Timepoint)
colData$Factor2<-paste(colData$Treatment, colData$Timepoint)


TPM<-as.data.frame(counts)
genes<-row.names(TPM)
TPM<-gather(TPM, key="Sample", value="TPM")
TPM$GeneID<-genes
TPM<-merge(TPM, colData, by="Sample")
TPM$Factor<-paste(TPM$GeneID, TPM$Rep, TPM$Timepoint, TPM$Genotype)
TPM2<-TPM[,c(9, 2, 5)]
TPM3<-spread(TPM2, key="Treatment", value="TPM")
TPM3$FoldChange<-TPM3$T/TPM3$C
metadata<-TPM[TPM$Treatment=="T",]
metadata<-metadata[,c(9, 3, 4, 6, 7)]
TPM3<-merge(TPM3, metadata, by="Factor")

all_significant$Factor<-paste(all_significant$Timepoint, all_significant$Regulation)
test<-list()
for(factor in unique(all_significant$Factor)){
  data<-all_significant[(all_significant$Factor==factor),]
  time<-unique(data$Timepoint)
  DEGs<-as.character(unique(data$gene))
  degs_tpm<-subset(TPM3, TPM3$GeneID %in% DEGs)
  degs_tpm<-subset(degs_tpm, degs_tpm$Timepoint %in% time)
  list<-list()
  for(gene in unique(degs_tpm$GeneID)){
    data<-degs_tpm[(degs_tpm$GeneID==gene),]
    kt<-kruskal.test(FoldChange ~ Genotype, data)[[3]]
    kt<-as.data.frame(cbind(unique(data$GeneID), kt))
    list[[(length(list)+1)]]<-kt
  }
  krusk<-as.data.frame(do.call("rbind", list))
  krusk$data<-paste(factor)
  test[[(length(test)+1)]]<-krusk
}

all<-as.data.frame(do.call("rbind", test))
write.csv(all, file="~/Documents/bmc/Data/compare_genotypes.csv")

