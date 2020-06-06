library(VennDiagram)

all_significant <- 
  read.csv("~/Documents/S L Fronteirs/all_significant_parents.csv")
all_significant<-
  all_significant[(all_significant$Cultivar=="Stigg" | all_significant$Cultivar=="Longbow"),]
all_significant<- all_significant[- grep("LC", all_significant$ID),]

# ==================================================================
Shared_genes_for_venn <- read.csv("~/Documents/bmc/Data/DE_genes/Shared/Shared_genes_for_venn.csv")



# subset genes that are DE in Stigg and in Longbow
Stigg<-all_significant[(all_significant$Cultivar=="Stigg"),]
SU<-Stigg[Stigg$Regulation=="Up",]
Longbow<-all_significant[(all_significant$Cultivar=="Longbow"),]
LU<-Stigg[Longbow$Regulation=="Up",]

Stigg<-as.character(unique(Stigg$ID))
Longbow<-as.character(unique(Longbow$ID))

shared<-subset(Longbow, Longbow %in% Stigg)
shared2<-subset(Stigg, Stigg %in% Longbow)
shared_DEG<-subset(all_significant, all_significant$ID %in% shared)

likes <- function(animals) {
  ppl <- Shared_genes_for_venn
  names(ppl) <- colnames(Shared_genes_for_venn)
  for (i in 1:length(animals)) {
    ppl <- subset(ppl, ppl[animals[i]] == T)
  }
  nrow(ppl)
}

plotAnimals <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(likes(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                            likes(a[2:3]), likes(a[c(1, 3)]), likes(a), euler.d = TRUE, ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                          likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                          likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2, 
                                                                                     4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), euler.d = TRUE, likes(a), ...)
  }
  if (!exists(out)) 
    out <- Oops
  return(out)
}

plotAnimals(c("LongbowUp", "LongbowDown", "StiggUp", "StiggDown"),
            category = c("LongbowUp", "LongbowDown", "StiggUp", "StiggDown"), 
            alpha=0.5, scaled=F, fill=c("lightgoldenrod1", "darkolivegreen1", "wheat3", "navajowhite2"), 
            cat.fontfamily =rep("sans", 4),  cex = rep(1.5, 15), 
            fontfamily = rep("sans", 15), cat.cex=rep(1.5, 4))


