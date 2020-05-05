library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
allowWGCNAThreads()

# Import kallisto files with txi 
samples<-as.data.frame(dir("Samples/"))
files <- file.path(samples$`dir("Samples/")`, "abundance.h5")
names(files) <- paste0(samples)
all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-colData %>%
  slice(match(order, Sample))
# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Treatment + Timepoint + Genotype + Treatment:Genotype)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Timepoint", "Genotype", "Rep"), returnData=TRUE)
vst_counts<-assay(vsd)


