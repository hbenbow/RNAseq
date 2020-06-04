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
rm(order)


# check that order of samples in metadata and txi object are the same
order<-colData$Sample
colData<-colData %>%
  slice(match(order, Sample))

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, 
                                colData, ~ Treatment + 
                                  Timepoint + Genotype + Treatment:Genotype)
vst<-varianceStabilizingTransformation(dds)
vds<-assay(vst)


# ==================================================================================
# Stigg specific
Stigg<-all_significant[(all_significant$Cultivar=="Stigg"),]
Stigg<-as.character(Stigg$ID)
Longbow<-all_significant[(all_significant$Cultivar=="Longbow"),]
Longbow<-as.character(Longbow$ID)

stigg_u<-subset(Stigg, !(Stigg %in% Longbow))
stigg_u<-as.character(unique(stigg_u))
stigg_u_exp<-subset(all_significant, all_significant$ID %in% stigg_u)
tpm<-txi.kallisto.tsv$counts
stigg_u_TPM<-vds[stigg_u,]
df<-stigg_u_TPM
df<-scale(stigg_u_TPM)

# ==================================================================================
# elbow method
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:5

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# ==================================================================================
# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(df, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")
fviz_nbclust(df, kmeans, method = "silhouette")
# ==================================================================================

# compute gap statistic
set.seed(123)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Print the result
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

# ==================================================================================
# choose value of k vased on three above methods
final <- kmeans(df, 3, nstart = 25)

fviz_cluster(final, data = df)

# ==================================================================================
DEGs<-as.character(unique(all_significant$gene))
deg_TPM<-vds[DEGs,]
df<-deg_TPM


pc<-prcomp(df)
plot(pc$x)

df<-scale(deg_TPM)
d <- dist(df, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, hang = 0, cex=.4)

# name and shame the outliers
groupindexes <- cutree(hc1, h = 5) # cut the tree at 5
table(groupindexes)
# this leaves 565 genes in one main cluster
df<-df[groupindexes == 1,] # make a new matrix using only genes from group 1

d <- dist(df, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, hang = -1, cex=.4)



