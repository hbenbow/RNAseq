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
library("genefilter")
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
rm(order)

hc<-as.character(all_expressed_genes$GeneID)
# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~ Treatment + 
                                  Timepoint + Genotype + Treatment:Genotype)
vst<-varianceStabilizingTransformation(dds)
vds<-assay(vst)
str(assay(vst))
hc<-vst[hc,]

topVarGenes <- head(order(rowVars(assay(hc)), decreasing = TRUE), 50)
mat  <- assay(vst)[ shared_160, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(hc)[, c("Rep", "Treatment", "Genotype", "Timepoint")])
pheatmap(mat, annotation_col = anno)

{shared_160<-c("TraesCS1A02G011500.1",
               "TraesCS1A02G234000.1",
               "TraesCS1A02G236300.1",
               "TraesCS1A02G358400.1",
               "TraesCS1A02G370600.1",
               "TraesCS1A02G370700.1",
               "TraesCS1A02G376700.1",
               "TraesCS1A02G376900.1",
               "TraesCS1A02G385400.3",
               "TraesCS1B02G129600.3",
               "TraesCS1B02G298900.1",
               "TraesCS1B02G337800.1",
               "TraesCS1D02G053700.2",
               "TraesCS1D02G085600.1",
               "TraesCS1D02G145500.4",
               "TraesCS1D02G159400.8",
               "TraesCS1D02G191600.1",
               "TraesCS1D02G288300.1",
               "TraesCS1D02G328400.1",
               "TraesCS1D02G376800.1",
               "TraesCS1D02G383900.1",
               "TraesCS1D02G389200.1",
               "TraesCS2A02G166600.2",
               "TraesCS2A02G227100.2",
               "TraesCS2A02G283600.1",
               "TraesCS2A02G316400.1",
               "TraesCS2A02G416900.1",
               "TraesCS2A02G420900.1",
               "TraesCS2A02G421800.1",
               "TraesCS2A02G460500.1",
               "TraesCS2B02G232300.1",
               "TraesCS2B02G236500.1",
               "TraesCS2B02G253400.1",
               "TraesCS2B02G369200.1",
               "TraesCS2B02G495700.2",
               "TraesCS2B02G593000.2",
               "TraesCS2D02G059200.1",
               "TraesCS2D02G279500.1",
               "TraesCS2D02G315600.2",
               "TraesCS2D02G321200.1",
               "TraesCS2D02G346900.2",
               "TraesCS2D02G386200.2",
               "TraesCS2D02G436200.1",
               "TraesCS2D02G438400.1",
               "TraesCS2D02G459200.1",
               "TraesCS2D02G485100.1",
               "TraesCS2D02G580700.1",
               "TraesCS2D02G598100.1",
               "TraesCS3A02G056800.1",
               "TraesCS3A02G067500.1",
               "TraesCS3A02G086700.4",
               "TraesCS3A02G087800.2",
               "TraesCS3A02G092000.2",
               "TraesCS3A02G395000.2",
               "TraesCS3A02G432800.1",
               "TraesCS3A02G492600.1",
               "TraesCS3A02G501200.1",
               "TraesCS3B02G014000.1",
               "TraesCS3B02G080100.1",
               "TraesCS3B02G126300.1",
               "TraesCS3B02G152300.3",
               "TraesCS3B02G238800.7",
               "TraesCS3B02G260400.1",
               "TraesCS3B02G383900.1",
               "TraesCS3B02G468800.1",
               "TraesCS3B02G565500.1",
               "TraesCS3D02G010700.1",
               "TraesCS3D02G017700.2",
               "TraesCS3D02G034200.1",
               "TraesCS3D02G034200.2",
               "TraesCS3D02G094300.1",
               "TraesCS3D02G167200.1",
               "TraesCS3D02G177400.1",
               "TraesCS3D02G214100.2",
               "TraesCS3D02G311100.1",
               "TraesCS3D02G457200.1",
               "TraesCS3D02G500100.1",
               "TraesCS4A02G136500.2",
               "TraesCS4B02G064900.1",
               "TraesCS4B02G094300.1",
               "TraesCS4B02G114700.1",
               "TraesCS4B02G316800.1",
               "TraesCS4B02G376300.1",
               "TraesCS4D02G063500.1",
               "TraesCS4D02G091100.1",
               "TraesCS4D02G169500.2",
               "TraesCS4D02G273200.4",
               "TraesCS5A02G020000.1",
               "TraesCS5A02G048400.1",
               "TraesCS5A02G068200.1",
               "TraesCS5A02G070300.1",
               "TraesCS5A02G079200.1",
               "TraesCS5A02G129500.1",
               "TraesCS5A02G138800.1",
               "TraesCS5A02G171600.1",
               "TraesCS5A02G175200.2",
               "TraesCS5A02G179200.1",
               "TraesCS5A02G311795.1",
               "TraesCS5A02G377000.5",
               "TraesCS5A02G478100.1",
               "TraesCS5A02G555100.1",
               "TraesCS5B02G138000.1",
               "TraesCS5B02G262100.1",
               "TraesCS5B02G448200.1",
               "TraesCS5B02G491000.1",
               "TraesCS5D02G007700.1",
               "TraesCS5D02G009100.1",
               "TraesCS5D02G130200.2",
               "TraesCS5D02G252700.1",
               "TraesCS5D02G414000.1",
               "TraesCS5D02G457200.1",
               "TraesCS5D02G491700.1",
               "TraesCS6A02G143500.1",
               "TraesCS6A02G179700.1",
               "TraesCS6A02G273000.1",
               "TraesCS6A02G298100.1",
               "TraesCS6A02G306800.1",
               "TraesCS6A02G321100.1",
               "TraesCS6A02G359000.2",
               "TraesCS6A02G381600.1",
               "TraesCS6B02G093300.2",
               "TraesCS6B02G095100.1",
               "TraesCS6B02G314300.4",
               "TraesCS6B02G323000.1",
               "TraesCS6B02G361500.1",
               "TraesCS6B02G446700.1",
               "TraesCS6D02G240000.1",
               "TraesCS6D02G257100.1",
               "TraesCS6D02G264800.1",
               "TraesCS6D02G289800.1",
               "TraesCS7A02G090300.1",
               "TraesCS7A02G096600.2",
               "TraesCS7A02G127700.3",
               "TraesCS7A02G175200.2",
               "TraesCS7A02G201300.1",
               "TraesCS7A02G201500.1",
               "TraesCS7A02G232700.2",
               "TraesCS7B02G108300.1",
               "TraesCS7B02G113800.1",
               "TraesCS7B02G157800.2",
               "TraesCS7B02G160200.2",
               "TraesCS7B02G414500.1",
               "TraesCS7B02G425800.1",
               "TraesCS7B02G432600.1",
               "TraesCS7D02G068200.1",
               "TraesCS7D02G086600.1",
               "TraesCS7D02G106200.1",
               "TraesCS7D02G177700.1",
               "TraesCS7D02G204600.1",
               "TraesCS7D02G204700.1",
               "TraesCS7D02G382200.2",
               "TraesCS7D02G384300.1",
               "TraesCS7D02G394900.1",
               "TraesCS7D02G410000.1",
               "TraesCS7D02G472100.1",
               "TraesCS7D02G493000.1",
               "TraesCSU02G049500.1",
               "TraesCSU02G087100.1",
               "TraesCSU02G187200.1",
               "TraesCSU02G198000.1")}