### 评估基因的富集结果

步骤：
对scPATIM识别出的遗传关联基因（TRGs）进行GO/KEGG通路富集分析。
验证这些通路是否与肿瘤免疫微环境（如T细胞耗竭、抗原呈递）直接相关。
对比其他方法（如MAGMA、FUMA）的富集结果，证明scPATIM的特异性。

```R
#scPagwas
magma_file="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/DLBCL_gene_pvalue.csv"
scpagwas_file="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/scPagwas_gene_PCC.csv"
scpatim_file="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/scPatim_gene_PCC.csv"
##读取结果
library(WebGestaltR)
setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/")
gene_df<-read.csv(scpatim_file)

refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
sig_genes<-gene_df$X[gene_df$adj_pvalue<=0.05]
outputDirectory <- getwd()
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=sig_genes,interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory=outputDirectory, projectName="scPatime_genes")

results<-read.table("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/Project_scPatime_genes/enrichment_results_scPatime_genes.txt",header = T,sep = "\t")
#绘制结果，选择前10个通路柱状图，根据p值排序，绘制p值的-log10，颜色根据p值大小
library(ggplot2)
library(dplyr)
results<-results[order(results$FDR),]
results<-results[1:10,]
results$FDR<-as.numeric(results$FDR)
results$logFDR<- -log10(results$FDR)
results<-na.omit(results)
results$description <- factor(results$description, levels = results$description[order(results$logFDR, decreasing = FALSE)])
pdf("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/scPatime_genes_enrichment.pdf",width=5,height=3)
ggplot(results,aes(x=description,y=logFDR))+geom_bar(stat="identity", fill = "skyblue")+coord_flip()+theme_minimal()+labs(title="scPatime_genes",x="Term",y="-log10(logFDR)")
dev.off()

##magma
gene_df<-read.csv(magma_file)
sig_genes<-gene_df$gene_symbol[gene_df$P<=0.05]
outputDirectory <- getwd()
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=sig_genes,interestGeneType="genesymbol", referenceGeneFile=refFile,
fdrThr=0.5,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory=outputDirectory, projectName="magma_genes")

##scPagwas
gene_df<-read.csv(scpagwas_file)
sig_genes<-gene_df$X[gene_df$PCC>=0.1]
outputDirectory <- getwd()
WebGestaltR(enrichMethod="ORA", organism="hsapiens",
enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=sig_genes,interestGeneType="genesymbol", referenceGeneFile=refFile,
referenceGeneType="genesymbol", isOutput=TRUE,
outputDirectory=outputDirectory, projectName="scPagwas_genes")
#读取结果
results<-read.table("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/Project_scPagwas_genes/enrichment_results_scPagwas_genes.txt",header = T,sep = "\t")
#绘制结果，选择前10个通路柱状图，根据p值排序，绘制p值的-log10，颜色根据p值大小
results<-results[order(results$FDR),]
results<-results[1:10,]
results$FDR<-as.numeric(results$FDR)
results$logFDR<- -log10(results$FDR)
results<-na.omit(results)
results$description <- factor(results$description, levels = results$description[order(results$logFDR, decreasing = FALSE)])
pdf("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/scPagwas_genes_enrichment.pdf",width=5,height=3)
ggplot(results,aes(x=description,y=logFDR))+geom_bar(stat="identity", fill = "skyblue")+coord_flip()+theme_minimal()+labs(title="scPagwas_genes",x="Term",y="-log10(logFDR)")
dev.off()



```