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
counts<-counts1
colData$TxGxT<-paste(colData$Rep, colData$Genotype, colData$Timepoint)
colData$Factor2<-paste(colData$Treatment, colData$Timepoint)
for(regulation in c("Up", "Down")){
  for(time in c(6, 24, 48, 96)){
    for(g in c("Stigg", "Longbow")){
      data<-all_significant[(all_significant$Regulation==regulation),]
      data<-data[(data$Timepoint==time),]
      data<-data[(data$Cultivar==g),]
      DEGs<-as.character(unique(data$gene))
      deg_TPM<-counts[DEGs,]
      deg_TPM<-as.data.frame(deg_TPM)
      genes<-row.names(deg_TPM)
      deg_TPM<-gather(deg_TPM, key="Sample", value="TPM")
      deg_TPM$GeneID<-genes
      deg_TPM<-merge(deg_TPM, colData, by="Sample")
      deg_TPM<-deg_TPM[(deg_TPM$Timepoint==time),]
      colData2<-colData[(colData$Treatment=="T"),]
      colData2<-colData2[,c(8, 5, 4, 2, 9)]

        data<-deg_TPM
        data$factor<-paste(data$GeneID,data$Rep )
        data<-data[,c(12,1,2)]
        data<-spread(data, key="Treatment", value="TPM")
        data$difference<-as.numeric(data$T/data$C)
        data$TxGxT<-paste(i)
        list[[(length(list)+1)]]<-data
      }
      normalised<-as.data.frame(do.call("rbind", list))
      normalised<-normalised[,c(1, 4, 5)]
      normalised<-left_join(normalised, colData2)
    }
  }
}
normalised$factor<-paste(normalised$GeneID, normalised$Genotype)
normalised<-na.omit(normalised)
normalised<-normalised[!(normalised$difference=="Inf"),]
res.aov3 <- aov(difference ~ factor,
                data = normalised)

tuk<-TukeyHSD(res.aov3, which = "factor")



}



ks.test(normalised$difference, "pnorm", 
        alternative = "two.sided")


for(i in unique(normalised$Timepoint)){
  data<-normalised[normalised$Timepoint==i,]
  res.aov3 <- aov(difference ~ factor,
                  data = data)
  tuk<-TukeyHSD(res.aov3, which = "factor")
  tukey<-as.data.frame(tuk$factor)
  sig<-tukey[tukey$padj < 0.05,]
  assign(paste("tuk",i), sig)
}

res.aov2 <- aov(difference ~ factor + Genotype + factor:Genotype, 
                data = normalised)





heatmap(df)

sampleDists <- dist(t(df))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(df), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)





# ==================================================================================
# make matrix if folc change values of each gene
all_significant$Factor<-paste(all_significant$Cultivar, all_significant$Timepoint)
exp_mat<-all_significant[,c(9, 20, 2)]
exp_mat<-spread(exp_mat, key="Factor", value="log2FoldChange", fill=0)
row.names(exp_mat)<-exp_mat$ID
exp_mat$ID<-NULL
exp_mat<-as.matrix(exp_mat)
df<-scale(exp_mat)
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

deg_expressed<-as.character(unique(all_significant$ID))
deg_expressed<-vds[deg_expressed,]
df<-scale(deg_expressed)

# ==================================================================================
# shared genes
all_significant <- 
  read.csv("~/Documents/S L Fronteirs/all_significant_parents.csv")
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
all_significant<- all_significant[- grep("LC", all_significant$ID),]

# ==================================================================
# subset genes that are DE in Stigg and in Longbow
Stigg<-all_significant[(all_significant$Cultivar=="Stigg"),]
Longbow<-all_significant[(all_significant$Cultivar=="Longbow"),]
Stigg<-as.character(unique(Stigg$ID))
Longbow<-as.character(unique(Longbow$ID))

shared<-subset(Longbow, Longbow %in% Stigg)
shared2<-subset(Stigg, Stigg %in% Longbow)
shared_DEG<-subset(all_significant, all_significant$ID %in% shared)
shared_exp<-vds[shared2,]
df<-scale(shared_exp)
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
                    K.max =50, B = 50)
# Print the result
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

# ==================================================================================
# choose value of k vased on three above methods
final <- kmeans(df, 2, nstart = 25)

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



