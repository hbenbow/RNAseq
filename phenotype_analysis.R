library(ggplot2)

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

pairwise.wilcox.test(Stigg_Longbow_phenotypes$Chlorosis,Stigg_Longbow_phenotypes$Factor,
                     p.adjust.method = "BH")
pairwise.wilcox.test(Stigg_Longbow_phenotypes$Pycnidia,Stigg_Longbow_phenotypes$Factor,
                     p.adjust.method = "BH")

# Make tables of mean, sd, se and n of phenotype data
pycnidia<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="mean")
pycnidia$SD<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="sd")[,3]
pycnidia$n<-aggregate(Stigg_Longbow_phenotypes$Pycnidia, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="length")[,3]
pycnidia$SE<-pycnidia$SD/sqrt(pycnidia$n)
colnames(pycnidia)<-c("Genotype", "Treatment", "Pycnidia", "SD", "n", "SE")

Chlorosis<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="mean")
Chlorosis$SD<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="sd")[,3]
Chlorosis$n<-aggregate(Stigg_Longbow_phenotypes$Chlorosis, list(Stigg_Longbow_phenotypes$Name, Stigg_Longbow_phenotypes$Treatment), FUN="length")[,3]
Chlorosis$SE<-Chlorosis$SD/sqrt(Chlorosis$n)
colnames(Chlorosis)<-c("Genotype", "Treatment", "Chlorosis", "SD", "n", "SE")

ggplot(Chlorosis, aes(x=Genotype, y=Chlorosis)) + 
  geom_bar(aes(fill=Treatment),position="dodge", stat="identity", alpha=0.7) + 
  theme_classic() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=Chlorosis - SE, ymax= Chlorosis+SE, group=Treatment),position=position_dodge(width=0.9), width=0.5) + 
  scale_fill_manual(values=c("black","grey50"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Percentage leaf area bearing chlorosis") + 
  theme(text = element_text(size=20, colour="black"), axis.text.x = element_text(colour="black")) 

ggplot(pycnidia, aes(x=Genotype, y=Pycnidia)) + 
  geom_bar(aes(fill=Treatment),position="dodge", stat="identity", alpha=0.7) + 
  theme_classic() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=Pycnidia - SE, ymax= Pycnidia+SE, group=Treatment),position=position_dodge(width=0.9), width=0.5) + 
  scale_fill_manual(values=c("black","grey50"), labels=c("Tween20", expression(paste(italic("Z. tritici"))))) + 
  ylab("Percentage leaf area bearing pycnidia") + 
  theme(text = element_text(size=20, colour="black"), axis.text.x = element_text(colour="black")) 

