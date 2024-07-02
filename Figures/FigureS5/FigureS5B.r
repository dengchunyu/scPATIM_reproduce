cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
#'finngen_r7_c3_bladder_exallc',
gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')

library(ggplot2)
len_1<-c()
for(i in 1:14){
    cancer<-cancers[i]
    gwas<-gwass[i]
    PCC<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_gene_PCC.csv"))
    len_1[i]<-sum(PCC$weight_pcc>0)
}

dd<-data.frame(number1=unlist(len_1),cancer=cancers)
dd <- dd[order(-dd$number1), ]
#绘制柱状图，横坐标为癌症类型，纵坐标为number1，对每一个柱子的上方加入数字标签，并且在y=500处加入一条虚线
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result")
pdf("barplot_gene_number.pdf",width=8,height=6)
ggplot(dd, aes(x=reorder(cancer, -number1), y=number1)) +
  geom_bar(stat="identity", fill="#006769") +
  geom_text(aes(label=number1), vjust=-0.3, color="#151515", size=3.5) +
  geom_hline(yintercept=500, linetype="dashed", color="#E6FF94") +
  labs(x="Cancer Type", y="Numbers of gene(weighted_pcc>0)") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()