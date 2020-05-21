library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
allowWGCNAThreads()

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files
samples<-as.data.frame(dir("Samples/"))
files <- file.path(samples$`dir("Samples/")`, "abundance.h5")
names(files) <- paste0(samples)
all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
# ==================================================================================

# if already have a txi object, load it with the metadata (colData)
load("~/Documents/bmc/Data/txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype + Treatment:Genotype)


# Identify numbers of expressed transcripts
# dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized=TRUE) >= .5 ) >=2
# dds <- dds[idx,]
# This would filter out genes where there are fewer 
# than 6 samples with normalized counts greater than or equal to .5. 
# Chose 6 as it is 1/8 of the samples (i.e. just treated in one genotype).
vst<-varianceStabilizingTransformation(dds)

expressed_genes$GeneID<-row.names(expressed_genes)
expressed_genes<- expressed_genes[- grep("LC", expressed_genes$GeneID),]
expressed_genes<-expressed_genes[,c(49, 1:48)]
expressed_genes_long<-expressed_genes %>% gather(Sample, TPM, 2:49)
all_wheat_genes<-merge(expressed_genes_long, colData, by="Sample")

test<-spread(test, Rep, TPM)
test<-aggregate(all_wheat_genes$TPM, by=list(all_wheat_genes$TPM, all_wheat_genes$GeneID), FUN="sum")

expressed_genes_filtered<-expressed_genes_long[(expressed_genes_long$TPM>=0.5),]
length(unique(expressed_genes_filtered$GeneID))
