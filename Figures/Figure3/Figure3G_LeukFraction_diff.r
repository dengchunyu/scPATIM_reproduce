leuk.frac<-read.table("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime_final/Figures/Figure3/TCGA_all_leuk_estimate.masked.20170107.tsv")
library(ggplot2)
library(dplyr)
library("ggsignif")
#查看CancerType的类别
colnames(leuk.frac)<-c('CancerType','samples','LeukFrac')
unique(leuk.frac$CancerType)
#对unique(leuk.frac$CancerType)分类
ad<-list(
"PCPG"=c("ACC","PCPG"),
"BLCA"=c("BLCA") ,
"BRCA"=c("BRCA" ),
"CESC"=c( "CESC"),
"COREAD"=c("COAD", "READ"),
"ESCA" =c("ESCA"),
"BC"=c("GBM" , "LGG"),
"HNSC" =c("HNSC"),
"RCC"=c("KICH" ,"KIRC","KIRP" ),
"LIHC" =c("LIHC"),
"NCSLC"=c("LUAD", "LUSC"),  
"OV"  =c("OV" ),
"PAAD" =c("PAAD"),
"PRAD" =c("PRAD"),
"SKCM" =c("SKCM"),
"STAD" =c("STAD"),
"THCA" =c("THCA"),
"UCEC"=c("UCEC", "UCS"),
"EAC"=c("UVM"),
"TGCT"=c("TGCT"))
M1<-c("BRCA","LUAD", "LUSC","SKCM","STAD","THCA","UCEC")
M2<-c("COAD", "READ","ESCA","PAAD")
M3<-c("OV" ,"PRAD","HNSC","KICH" ,"KIRC","KIRP","LIHC")
leuk.frac<-leuk.frac[leuk.frac$CancerType %in% c(M1,M2,M3),]
leuk.frac$Type<-"Module3"
leuk.frac$Type[leuk.frac$CancerType %in% M1]<-"Module1"
leuk.frac$Type[leuk.frac$CancerType %in% M2]<-"Module2"
leuk.frac$LeukFrac<-as.numeric(leuk.frac$LeukFrac)
p <- ggplot(leuk.frac, aes(x = Type, y = LeukFrac, fill = Type)) +
  geom_violin() + 
  geom_boxplot(width = 0.1) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3) +
  theme_minimal()+
  scale_fill_manual(values=c("#DF826C","#9FBB73","#7C93C3"))+
  geom_signif(step_increase = c(0.1,0.1,0.1), comparisons=list(c("Module1","Module2"),c("Module1","Module3"),"Module2","Module3"),
      map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05), textsize=2,tip_length=0.03)

ggsave(p,file="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime_final/Figures/Figure3/Figure3G_LeukFrac_diff_module.pdf",width=6,height=5)

pairwise_comparisons <- pairwise.wilcox.test(leuk.frac$LeukFrac, leuk.frac$Type, p.adjust.method = "BH")

# 提取p值
p_values <- pairwise_comparisons$p.value
#p_values
#             Module1      Module2
#Module2 4.742762e-05           NA
#Module3 3.426634e-33 4.985454e-06