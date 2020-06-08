stigg_spec_PP<-as.character(c("TraesCS6A02G270800.10", 
                              "TraesCS3B02G238800.6", 
                              "TraesCS5A02G225200.1", 
                              "TraesCS5B02G433300.2", 
                              "TraesCS5D02G073800.1", 
                              "TraesCS5A02G186500.3", 
                              "TraesCS2A02G391100.1",
                              "TraesCS2D02G508800.3" ,
                              "TraesCS6B02G128000.1"))

stigg_spec_PP<-c("TraesCS1A02G001900.9",
"TraesCS4D02G209200.2",
"TraesCS4D02G147100.1",
"TraesCS3D02G309800.1",
"TraesCS2B02G585100.2",
"TraesCS3D02G309800.1",
"TraesCS3D02G309800.1",
"TraesCS2B02G585100.2",
"TraesCS7D02G495200.1",
"TraesCS7D02G495200.1",
"TraesCS4A02G013400.1")

development<-c("TraesCS2A02G420900.1", "TraesCS7B02G432600.1")
stigg_spec_PP<-development

other<-c("TraesCS2B02G585100.2",
"TraesCS3D02G309800.1",
"TraesCS4A02G013400.1",
"TraesCS4D02G209200.2",
"TraesCS7D02G495200.1")
stigg_spec_PP_other

mostactive<-c("TraesCS1D02G053700.2",
              "TraesCS3B02G238800.7",
              "TraesCS3D02G177400.1",
              "TraesCS5A02G068200.1",
              "TraesCS5A02G179200.1",
              "TraesCS5D02G130200.2",
              "TraesCS6A02G273000.1")
stigg_spec_PP<-mostactive

transport_stigg<-c("TraesCS2A02G391100.1",
"TraesCS2D02G508800.3" ,
"TraesCS6B02G128000.1")
stigg_spec_PP<-transport_stigg

transport<-c("TraesCS7A02G516600.2", "TraesCS7B02G311400.2", "TraesCS3A02G036600.2")
stigg_spec_PP<-transport

pp_stigg_gnees<-c("TraesCS3B02G238800.7", "TraesCS3B02G424200.3", "TraesCS5B02G528300.2")
stigg_spec_PP<-pp_stigg_gnees

pp_fc<-subset(all_significant, all_significant$ID %in% stigg_spec_PP)
pp_exp<-as.data.frame(vds[stigg_spec_PP,])
ppcor<-cor(t(pp_exp))
pp_exp<-gather(pp_exp, key=Sample, value="TPM")
pp_genes<-gather(pp_exp, key="Sample", value="TPM")
pp_genes$GeneID<-paste(stigg_spec_PP)
pp_genes<-merge(pp_genes, colData, by="Sample")
pp_genes2<-aggregate(pp_genes$TPM, by=list(pp_genes$GeneID, pp_genes$Treatment, pp_genes$Timepoint, pp_genes$Genotype), FUN="mean")
pp_genes2$sd<-aggregate(pp_genes$TPM, by=list(pp_genes$GeneID, pp_genes$Treatment, pp_genes$Timepoint, pp_genes$Genotype), FUN="sd")$x
pp_genes2$n<-aggregate(pp_genes$TPM, by=list(pp_genes$GeneID, pp_genes$Treatment, pp_genes$Timepoint, pp_genes$Genotype), FUN="length")$x
pp_genes2$SE<-pp_genes2$sd/sqrt(pp_genes2$n)

pp_genes2$Group.1 = factor(pp_genes2$Group.1, levels=c("TraesCS6A02G270800.10", 
                                                       "TraesCS3B02G238800.6", 
                                                       "TraesCS5A02G225200.1", 
                                                       "TraesCS5B02G433300.2", 
                                                       "TraesCS5D02G073800.1", 
                                                       "TraesCS5A02G186500.3", 
                                                       "TraesCS2A02G391100.1",
                                                       "TraesCS2D02G508800.3" ,
                                                       "TraesCS6B02G128000.1")) 


ggplot(pp_genes2, aes(x=as.factor(Group.3), y=x, group=Group.2))+
  geom_line(aes(colour=Group.2), size=1, alpha=0.7) +
  facet_grid(Group.4~Group.1) +
  geom_errorbar(aes(ymin=x-SE,ymax=x+SE, colour=Group.2), width=.2, size=1, alpha=0.7)+
  theme_bw() +
  scale_colour_manual(values=c("grey30", "orange"), 
                      labels=c("Tween20", expression(paste(italic("Z. tritici")))),
                      labs(colour="Treatment")) +
  xlab("Hours post inoculation") +
  theme(text = element_text(size=20, colour='black'), axis.text.x = element_text(colour="black")) +
  ylab("Normalised transcript count")

ggplot(pp_fc, aes(x=as.factor(Timepoint), y=log2FoldChange)) +
  geom_col(aes(fill=ID), position=position_dodge(preserve="single"))+
  theme_classic()

