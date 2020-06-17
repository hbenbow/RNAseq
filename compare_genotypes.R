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
library(bestNormalize)
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


order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))
rm(order)
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~ Factor)
dds2<-DESeq(dds)
REST6<-as.data.frame(results(dds2, contrast=c("Factor", "T6Longbow", "T6Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
REST24<-as.data.frame(results(dds2, contrast=c("Factor", "T24Longbow", "T24Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
REST48<-as.data.frame(results(dds2, contrast=c("Factor", "T48Longbow", "T48Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
REST96<-as.data.frame(results(dds2, contrast=c("Factor", "T96Longbow", "T96Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RESC6<-as.data.frame(results(dds2, contrast=c("Factor", "C6Longbow", "C6Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RESC24<-as.data.frame(results(dds2, contrast=c("Factor", "C24Longbow", "C24Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RESC48<-as.data.frame(results(dds2, contrast=c("Factor", "C48Longbow", "C48Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RESC96<-as.data.frame(results(dds2, contrast=c("Factor", "C96Longbow", "C96Stigg"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))

REST6$Treatment<-"Treated"
REST24$Treatment<-"Treated"
REST48$Treatment<-"Treated"
REST96$Treatment<-"Treated"
RESC6$Treatment<-"Control"
RESC24$Treatment<-"Control"
RESC48$Treatment<-"Control"
RESC96$Treatment<-"Control"

REST6$Timepoint<-6
REST24$Timepoint<-24
REST48$Timepoint<-48
REST96$Timepoint<-96
RESC6$Timepoint<-6
RESC24$Timepoint<-24
RESC48$Timepoint<-48
RESC96$Timepoint<-96

all_compare_geno<-rbind(REST6, REST24, REST48, REST96, RESC6, RESC24, RESC48, RESC96)
all_sig_genotype<-all_compare_geno[(all_compare_geno$padj<0.05),]
all_sig_genotype<-all_sig_genotype[abs(all_sig_genotype$log2FoldChange)>.5,]
all_sig_genotype<-na.omit(all_sig_genotype)

assay<-assay(dds2)
vst<-varianceStabilizingTransformation(dds2)
vds<-assay(vst)


counts1<-vds
# counts1<-vds
counts<-counts1
colData$TxGxT<-paste(colData$Rep, colData$Genotype, colData$Timepoint)
colData$Factor2<-paste(colData$Treatment, colData$Timepoint)


TPM<-as.data.frame(counts1)
genes<-row.names(TPM)
TPM<-gather(TPM, key="Sample", value="TPM")
TPM$GeneID<-genes
TPM<-merge(TPM, colData, by="Sample")
TPM$Factor<-paste(TPM$GeneID, TPM$Rep, TPM$Timepoint, TPM$Genotype)
TPM2<-TPM[,c(9, 2, 5)]
TPM3<-spread(TPM2, key="Treatment", value="TPM")
TPM3$FoldChange<-(TPM3$T)/(TPM3$C)
metadata<-TPM[TPM$Treatment=="T",]
metadata<-metadata[,c(9, 3, 4, 6, 7)]
TPM4<-merge(TPM3, metadata, by="Factor")
TPM4$Log2FoldChange<-sqrt(TPM4$FoldChange)


all_significant$Factor<-paste(all_significant$Timepoint, all_significant$Regulation)
test<-list()
for(factor in unique(all_significant$Factor)){
  data<-all_significant[(all_significant$Factor==factor),]
  time<-unique(data$Timepoint)
  DEGs<-as.character(unique(data$gene))
  degs_tpm<-subset(TPM4, TPM4$GeneID %in% DEGs)
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

write.csv(all, file="~/Documents/bmc/Data/compare_genotypes2.csv")

all_foldchange<-subset(TPM4, TPM4$GeneID %in% all$V1)

