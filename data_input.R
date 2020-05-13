library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(plyr)
allowWGCNAThreads()

# Import kallisto files with txi 
# ==================================================================================
# Use these steps to import kallisto files
samples<-dir("~/Documents/bmc/Data/samples/") # where the directory 'samples'
# contains the kallisto output directories - 1 per sample.
files <- file.path(samples, "abundance.h5")
names(files) <- paste0(samples)
all(file.exists(files)) # need to navigate into samples directory for this to work!!
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.csv("~/Documents/bmc/Data/colData.csv")
setwd("~/Documents/bmc/Data/")
colData$Timepoint<-as.factor(colData$Timepoint)
# ==================================================================================


# if already have a txi object, load it with the metadata (colData)
load("~/Documents/bmc/Data/txi_and_colData.RData")

# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype + Treatment:Genotype)

# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Timepoint", "Genotype", "Rep"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Genotype, shape=Treatment), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black")) +
  scale_colour_manual(values=c("black","grey50"))

vst_counts<-assay(vsd)


