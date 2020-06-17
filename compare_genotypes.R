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
write.csv(all_sig_genotype, file="~/Documents/bmc/Data/all_significant_by_genotype.csv")
# treated only DEGs
all_sig_genotype_genes<-as.character(unique(all_sig_genotype$row))
geno_and_treatment<-subset(all_significant, all_significant$ID %in% all_sig_genotype_genes)

all_significant$Factor<-paste(all_significant$Cultivar, all_significant$Timepoint, all_significant$Regulation)
Factor="Stigg 6 Down"

list<-list()
for(i in as.character(unique(all_significant$Factor))){
data<-all_significant[(all_significant$Factor==i),]
timepoint<-as.character(unique(data$Timepoint))
regulation<-as.character(unique(data$Regulation))
genos<-all_sig_genotype[(all_sig_genotype$Timepoint==timepoint),]
genos<-subset(genos, genos$row %in% data$ID)
chroms<-subset(all_significant, all_significant$ID %in% genos$row)
chroms<-chroms[!(duplicated(chroms$ID)),]
chroms<-chroms[,c(9, 12, 13, 14)]
genos<-merge(genos, chroms, by.x="row", by.y="ID", all.x=T, all.y=F)
t<-as.data.frame(table(genos$row, genos$Treatment, genos$Genome))
t<-spread(t, key="Var2", value="Freq")
tonly<-t[(t$Treated==1),]
tonly<-tonly[(tonly$Control==0),]
tonly<-spread(as.data.frame(table(tonly$Var3)), key="Var1", value="Freq")
tonly$Factor<-i
tonly$Treatment<-"Treated"
conly<-t[(t$Control==1),]
conly<-conly[(conly$Treated==0),]
conly<-spread(as.data.frame(table(conly$Var3)), key="Var1", value="Freq")
conly$Factor<-i
conly$Treatment<-"Control"
both<-t[(t$Control==1),]
both<-both[(both$Treated==1),]
both<-spread(as.data.frame(table(both$Var3)), key="Var1", value="Freq")
both$Factor<-i
both$Treatment<-"Both"
d2<-as.data.frame(rbind(tonly, conly, both))
list[[(length(list)+1)]]<-d2
}

table<-as.data.frame(do.call(rbind.data.frame, list))
write.csv(table, file="~/Documents/bmc/Data/DE_genes/table_of_degs.csv")

Stigg_DEGs<-as.character(unique(all_significant[(all_significant$Cultivar=="Stigg"),9]))
stigg_DEGs_geno<-subset(all_sig_genotype, all_sig_genotype$row %in% Stigg_DEGs)
chroms<-subset(all_significant, all_significant$ID %in% Stigg_DEGs)
chroms<-chroms[!(duplicated(chroms$ID)),]
chroms<-chroms[,c(9, 12, 13, 14)]
genos<-merge(stigg_DEGs_geno, chroms, by.x="row", by.y="ID", all.x=T, all.y=F)
