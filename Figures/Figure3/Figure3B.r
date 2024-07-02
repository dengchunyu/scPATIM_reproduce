library(stringr)
library(dplyr)
#library(GGally)
library(corrplot)
#'HeadandNeck',
cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass<-c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')

setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/")
#获取所有癌症得到的相关性基因
top_genelist<-list()
for(i in 1:15){
    cancer<-cancers[i]
    gwas<-gwass[i]
print(cancer)
PCC<-read.csv(paste0(cancer,"/",gwas,"_gene_PCC.csv"))
top_genelist[[i]]<- PCC$X[order(PCC$weight_pcc,decreasing=T)[1:500]]
}
#c("BladderCancer","GastricCancer","LiverCancer","ColorectalCancer")
top_genelist<-as.data.frame(top_genelist)
colnames(top_genelist)<-cancers

setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
# 假设top_gene_mat是一个数据框，每一列代表一个癌症包含的基因
# 创建一个空的矩阵来存储基因交集个数
n_cancers <- ncol(top_genelist)
intersection_matrix <- matrix(0, nrow = n_cancers, ncol = n_cancers)

# 计算每对癌症之间的基因交集个数
for (i in 1:(n_cancers - 1)) {
  for (j in (i + 1):n_cancers) {
    set1 <- top_genelist[, i]
    set2 <- top_genelist[, j]
    intersection_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    # 计算Jaccard系数
    jaccard_coefficient <- intersection_size / union_size
    intersection_matrix[i, j] <- jaccard_coefficient
    intersection_matrix[j, i] <- jaccard_coefficient
  }
}
colnames(intersection_matrix)<-cancers
rownames(intersection_matrix)<-cancers

load("D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/intersection_matrix.RData")

#install.packages("gplots")  # 如果未安装的话
library(gplots)
#setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
#corder<-c('Melanoma','OvarianCancer','BladderCancer','BreastCancer','LungCancer','Pancreatic','GastricCancer','HeadandNeck','EndometrialCarcinoma','ColorectalCancer','LiverCancer','EsophagealCancer','KidneyCancer','ProstateCancer','ThyroidCancer')
#intersection_matrix<-intersection_matrix[corder,corder]
# 绘制热图
col_fun<-colorRamp2(c(0,0.3), c("#F0F0F0","#3C486B"))

pdf("Figure3B_intersection_cancers_gene.pdf")
Heatmap(intersection_matrix, col = col_fun)
dev.off()