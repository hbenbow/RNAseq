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

setwd("~/Documents/bmc/Data/DE_genes/Shared/")

for(i in dir(pattern=".txt")){
  data<-read.delim(i)
  data$go_cat<-row.names(data)
  ggplot(data, aes(x=reorder(go_cat, GO), y=GO)) +
    geom_col() +
    coord_flip() +
    theme_classic()+
    theme(text = element_text(size=30, colour='black'), axis.text.x = element_text(colour="black"))+
    scale_x_discrete(labels=wrap_format(40))+
    ylab("Number of differentially expressed genes") +
    xlab("Biological process") +
    scale_y_continuous(labels = scales::number_format(accuracy = 1))
  ggsave(paste("~/Documents/bmc/Graphs/", i, ".pdf"))
}

LuSu<-Shared_genes_for_venn[(Shared_genes_for_venn$LongbowUp==1 & 
                               Shared_genes_for_venn$StiggUp==1 &
                               Shared_genes_for_venn$LongbowDown==0 &
                               Shared_genes_for_venn$StiggDown ==0),]
LuSu$Category<-"LuSu"

LuSd<-Shared_genes_for_venn[(Shared_genes_for_venn$LongbowUp==1 & 
                               Shared_genes_for_venn$StiggUp==0 &
                               Shared_genes_for_venn$LongbowDown==0 &
                               Shared_genes_for_venn$StiggDown ==1),]
LuSd$Category<-"LuSd"

LdSd<-Shared_genes_for_venn[(Shared_genes_for_venn$LongbowUp==0 & 
                               Shared_genes_for_venn$StiggUp==0 &
                               Shared_genes_for_venn$LongbowDown==1 &
                               Shared_genes_for_venn$StiggDown ==1),]
LdSd$Category<-"LdSd"

LdSu<-Shared_genes_for_venn[(Shared_genes_for_venn$LongbowUp==0 & 
                               Shared_genes_for_venn$StiggUp==1 &
                               Shared_genes_for_venn$LongbowDown==1 &
                               Shared_genes_for_venn$StiggDown ==0),]
LdSu$Category<-"Ldsu"

all<-rbind(LuSu, LuSd, LdSd, LdSu)

shared_exp<-subset(all_significant, all_significant$ID %in% all$Row.Labels)
shared_exp<-merge(shared_exp, all, by.x="ID", by.y="Row.Labels")

shared_exp<-shared_exp[,c(1, 3, 4, 7, 9, 10, 20)]





