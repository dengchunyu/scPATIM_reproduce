library(ggplot2)
library(dplyr)
library("ggsignif")
#读取TMB_OS.csv数据
setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime_final/Figures/Figure3/")
TMB_OS<-read.csv("TMB_OS.csv",header = T,sep = ",")
#处理TMB.median.interquartile中类似151(85-269)的数据，只保留151
TMB_OS$TMB.median<-gsub("\\(.*\\)","",TMB_OS$TMB.median.interquartile)
#将TMB.median.interquartile转换为数值型
TMB_OS$TMB.median<-as.numeric(TMB_OS$TMB.median)

ca<-intersect(c(M1,M2,M3),TMB_OS$tumor)
TMB_OS<-TMB_OS[TMB_OS$tumor %in% ca,]
M1<-c("BRCA","LUAD", "LUSC","SKCM","STAD","THCA","UCEC")
M2<-c("COAD", "READ","ESCA","PAAD")
M3<-c("OV" ,"PRAD","HNSC","KICH" ,"KIRC","KIRP","LIHC")
dd<-TMB_OS[,c("tumor","TMB.median")]
dd$Type<-"TRC-T"
dd$Type[dd$tumor %in% M2]<-"TRC-B"
dd$Type[dd$tumor %in% M3]<-"TRC-Myeloid"

dd<-dd[order(dd$TMB.median,decreasing=T),]
dd$tumor<-factor(dd$tumor,levels=dd$tumor)

p_tmb<-ggplot(data=dd, mapping=aes(x = tumor, y =TMB.median,fill=Type))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#DF826C","#9FBB73","#7C93C3"))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),legend.position = "top")+
  labs(y = "TMB")

ggsave(p_tmb,file="Figure3I_TMB_diff_module.pdf",width=6,height=5)


dNdS_escape<-read.csv("dNdS_escape.csv",header = T,sep = ",")
dNdS_escape<-dNdS_escape[,c("Tumor","Median.idN.dS.escaped","Median.IdN.dS.edited")]

ca<-intersect(c(M1,M2,M3),dNdS_escape$Tumor)
dNdS_escape<-dNdS_escape[dNdS_escape$Tumor %in% ca,]

dd<-dNdS_escape[,c("Tumor","Median.idN.dS.escaped")]
dd$Type<-"TRC-T"
dd$Type[dd$tumor %in% M2]<-"TRC-B"
dd$Type[dd$tumor %in% M3]<-"TRC-Myeloid"

dd<-dd[order(dd$Median.idN.dS.escaped,decreasing=T),]
dd$Tumor<-factor(dd$Tumor,levels=dd$Tumor)

p_dn<-ggplot(data=dd, mapping=aes(x = Tumor, y =Median.idN.dS.escaped,fill=Type))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#DF826C","#9FBB73","#7C93C3"))+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),legend.position = "top")+
  labs(y = "Median.idN.dS.escaped")

ggsave(p_dn,file="Figure3J_Median.idN.dS.escaped_diff_module.pdf",width=6,height=5)