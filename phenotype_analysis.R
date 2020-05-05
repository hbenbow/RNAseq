library(ggplot2)
library(multcompView)
        
# load phenotype data for Stigg and Longbow
Stigg_Longbow_phenotypes <- 
  read.csv("~/Documents/bmc/Data/Stigg_Longbow_phenotypes.csv")

# Create a new column of "factor" which describes genotype and treatment
Stigg_Longbow_phenotypes$Factor<-
  paste(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment)
# convert to factor
Stigg_Longbow_phenotypes$Factor<-as.factor(Stigg_Longbow_phenotypes$Factor)
# remove NAs
Stigg_Longbow_phenotypes<-na.omit(Stigg_Longbow_phenotypes)

# test for normality
ks.test(Stigg_Longbow_phenotypes$Pycnidia, "pnorm", 
        alternative = "two.sided")
shapiro.test(Stigg_Longbow_phenotypes$Pycnidia)

# define function that will calculate letters for pairwise comparisons
tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}

pp_chlor<-pairwise.wilcox.test(Stigg_Longbow_phenotypes$Chlorosis,Stigg_Longbow_phenotypes$Factor,
                     p.adjust.method = "BH")
pp_pyc<-pairwise.wilcox.test(Stigg_Longbow_phenotypes$Pycnidia,Stigg_Longbow_phenotypes$Factor,
                     p.adjust.method = "BH")

mymat_chlor<-tri.to.squ(pp_chlor$p.value)
mymat_pyc<-tri.to.squ(pp_pyc$p.value)

myletters_chlor<-as.data.frame(multcompLetters(mymat_chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
colnames(myletters_chlor)<-"Letters"
myletters_chlor$Factor<-row.names(myletters_chlor)
myletters_pyc<-as.data.frame(multcompLetters(mymat_pyc,compare="<=",threshold=0.05,Letters=letters)$Letters)
colnames(myletters_pyc)<-"Letters"
myletters_pyc$Factor<-row.names(myletters_pyc)

# Make tables of mean, sd, se and n of phenotype data
Chlorosis<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="mean")
Chlorosis$SD<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="sd")[,3]
Chlorosis$n<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="length")[,3]
Chlorosis$SE<-Chlorosis$SD/sqrt(Chlorosis$n)
colnames(Chlorosis)<-c("Genotype", "Treatment", "Chlorosis", "SD", "n", "SE")
Chlorosis$Factor<-paste(Chlorosis$Genotype, Chlorosis$Treatment)
Chlorosis<-merge(Chlorosis, myletters_chlor, by="Factor")


pycnidia<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="mean")
pycnidia$SD<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="sd")[,3]
pycnidia$n<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="length")[,3]
pycnidia$SE<-pycnidia$SD/sqrt(pycnidia$n)
colnames(pycnidia)<-c("Genotype", "Treatment", "Pycnidia", "SD", "n", "SE")
pycnidia$Factor<-paste(pycnidia$Genotype, pycnidia$Treatment)
pycnidia<-merge(pycnidia, myletters_pyc, by="Factor")



ggplot(Chlorosis, aes(x=Genotype, y=Chlorosis, group=Treatment)) + 
  geom_bar(aes(fill=Treatment),position="dodge", stat="identity", alpha=0.7) + 
  theme_classic() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=Chlorosis - SE, ymax= Chlorosis+SE, group=Treatment),position=position_dodge(width=0.9), width=0.5) + 
  scale_fill_manual(values=c("black","grey50"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Percentage leaf area bearing chlorosis") + 
  theme(text = element_text(size=20, colour="black"), axis.text.x = element_text(colour="black")) +
  geom_text(aes(x=Genotype, y=Chlorosis+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=5)+
  ylim(0,70)


ggplot(pycnidia, aes(x=Genotype, y=Pycnidia, group=Treatment)) + 
  geom_bar(aes(fill=Treatment),position="dodge", stat="identity", alpha=0.7) + 
  theme_classic() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=Pycnidia - SE, ymax= Pycnidia+SE, group=Treatment),position=position_dodge(width=0.9), width=0.5) + 
  scale_fill_manual(values=c("black","grey50"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Percentage leaf area bearing pycnidia") + 
  theme(text = element_text(size=20, colour="black"), axis.text.x = element_text(colour="black")) +
  geom_text(aes(x=Genotype, y=Pycnidia+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=5) +
  ylim(0,70)

save.image("~/Documents/bmc/Data/phenotype_data.RData")
