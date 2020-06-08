
colData_Longbow$L31<-MEs_Longbow$MEpaleturquoise
colData_Longbow$L17<-MEs_Longbow$MEgrey60
colData_Stigg$S20<-MEs_Stigg$MEroyalblue

ggplot(colData_Longbow, aes(x=Timepoint, y=L17)) + 
  geom_boxplot(aes(fill=Treatment)) + 
  theme_classic() +
  theme(text = element_text(size=40, colour='black'), axis.text.x = element_text(colour="black")) +
  xlab("Hours post inoculation") +
  ylab("Eigengene of Longbow module 17") +
  scale_fill_manual(values=c("grey30", "grey70"), labels=c("Control", "Treated")) +
  ylim(-.6,.6) +
  geom_hline(aes(yintercept=0))
ggsave("~/Documents/bmc/Graphs/MEL17.pdf")

ggplot(colData_Longbow, aes(x=Timepoint, y=L31)) + 
  geom_boxplot(aes(fill=Treatment)) + 
  theme_classic() +
  theme(text = element_text(size=40, colour='black'), axis.text.x = element_text(colour="black")) +
  xlab("Hours post inoculation") +
  ylab("Eigengene of Longbow module 31") +
  scale_fill_manual(values=c("grey30", "grey70"), labels=c("Control", "Treated")) +
  ylim(-.6,.6) +
  geom_hline(aes(yintercept=0))
ggsave("~/Documents/bmc/Graphs/MEL31.pdf")


ggplot(colData_Stigg, aes(x=Timepoint, y=S20)) + 
  geom_boxplot(aes(fill=Treatment)) + 
  theme_classic() +
  theme(text = element_text(size=40, colour='black'), axis.text.x = element_text(colour="black")) +
  xlab("Hours post inoculation") +
  ylab("Eigengene of Stigg module 20") +
  scale_fill_manual(values=c("grey30", "grey70"), labels=c("Control", "Treated")) +
  ylim(-.6,.6) +
  geom_hline(aes(yintercept=0))
ggsave("~/Documents/bmc/Graphs/MES20.pdf")

Longbow
MEpaleturquoise = L31
MEgrey60 = 17

MEroyalblue =S20

