
library(WebGestaltR)
setwd("D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
load("cancers_top500.RData")
refFile <- system.file("extdata", "referenceGenes.txt", package="WebGestaltR")
outputDirectory <- "D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/gene_enrichment"

for(i in colnames(gene_list)){
    pcc_gene<-gene_list[,i]
 WebGestaltR(enrichMethod="ORA", organism="hsapiens",
    enrichDatabase="pathway_Reactome", interestGene=pcc_gene,
    interestGeneType="genesymbol", referenceGeneFile=refFile,
    referenceGeneType="genesymbol", isOutput=TRUE,
    outputDirectory=outputDirectory, projectName=i,
    minNum = 5,
  maxNum = 100,
  sigMethod = "fdr",
  fdrMethod = "BH",
  fdrThr = 1,
  topThr = 100)
}

for(i in colnames(gene_list)){
    pcc_gene<-gene_list[,i]
 WebGestaltR(enrichMethod="ORA", organism="hsapiens",
    enrichDatabase="geneontology_Biological_Process_noRedundant", interestGene=pcc_gene,
    interestGeneType="genesymbol", referenceGeneFile=refFile,
    referenceGeneType="genesymbol", isOutput=TRUE,
    outputDirectory=outputDirectory, projectName=paste0(i,"_bp"),
    minNum = 5,
  maxNum = 100,
  sigMethod = "fdr",
  fdrMethod = "BH",
  fdrThr = 1,
  topThr = 100)
}
for(i in colnames(gene_list)){
    pcc_gene<-gene_list[,i]
 WebGestaltR(enrichMethod="ORA", organism="hsapiens",
    enrichDatabase="geneontology_Cellular_Component", interestGene=pcc_gene,
    interestGeneType="genesymbol", referenceGeneFile=refFile,
    referenceGeneType="genesymbol", isOutput=TRUE,
    outputDirectory=outputDirectory, projectName=paste0(i,"_cc"),
    minNum = 5,
  maxNum = 100,
  sigMethod = "fdr",
  fdrMethod = "BH",
  fdrThr = 1,
  topThr = 100)
}
for(i in colnames(gene_list)){
    pcc_gene<-gene_list[,i]
 WebGestaltR(enrichMethod="ORA", organism="hsapiens",
    enrichDatabase="geneontology_Molecular_Function", interestGene=pcc_gene,
    interestGeneType="genesymbol", referenceGeneFile=refFile,
    referenceGeneType="genesymbol", isOutput=TRUE,
    outputDirectory=outputDirectory, projectName=paste0(i,"_mf"),
    minNum = 5,
  maxNum = 100,
  sigMethod = "fdr",
  fdrMethod = "BH",
  fdrThr = 1,
  topThr = 100)
}

#geneontology_Cellular_Component
#geneontology_Molecular_Function
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(ggsci)
library(stringr)
outputDirectory <- "D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/gene_enrichment"

enrich_result_list<-list()
for(i in colnames(gene_list)){
  print(i)
    #判断是否有Project_+i文件夹，如果没有，那么跳过
    if(!file.exists(paste0(outputDirectory,"/Project_",i,"_bp/enrichment_results_",i,"_bp.txt"))){
      enrich_df_bp<-NULL
    }else{
    #读取Project_+i的文件夹中，enrichment_results_+i+.txt的文件
    enrich_df<-read.delim(paste0(outputDirectory,"/Project_",i,"_bp/enrichment_results_",i,"_bp.txt"),header = T)
    enrich_df<-enrich_df[enrich_df$FDR<0.05,]
    enrich_df$cancer<-i
    enrich_df$method<-"pos"
    enrich_df_bp<-enrich_df
    #enrich_result_list[[i]]<-enrich_df
    }
    if(!file.exists(paste0(outputDirectory,"/Project_",i,"_cc/enrichment_results_",i,"_cc.txt"))){
      enrich_df_cc<-NULL
    }else{
    #读取Project_+i的文件夹中，enrichment_results_+i+.txt的文件
    enrich_df<-read.delim(paste0(outputDirectory,"/Project_",i,"_cc/enrichment_results_",i,"_cc.txt"),header = T)
    enrich_df<-enrich_df[enrich_df$FDR<0.05,]
    enrich_df$cancer<-i
    enrich_df$method<-"CC"
    enrich_df_cc<-enrich_df
    #enrich_result_list[[i]]<-enrich_df
    }
    if(!file.exists(paste0(outputDirectory,"/Project_",i,"_mf/enrichment_results_",i,"_mf.txt"))){
      enrich_df_mf<-NULL
    }else{
    #读取Project_+i的文件夹中，enrichment_results_+i+.txt的文件
    enrich_df<-read.delim(paste0(outputDirectory,"/Project_",i,"_mf/enrichment_results_",i,"_mf.txt"),header = T)
    #enrich_df<-enrich_df[enrich_df$FDR<0.05,]
    enrich_df$cancer<-i
    enrich_df$method<-"MF"
    enrich_df_mf<-enrich_df
    #enrich_result_list[[i]]<-enrich_df
    }
    enrich_df<-rbind(enrich_df_bp,enrich_df_cc,enrich_df_mf)
    enrich_df<-enrich_df[order(enrich_df$FDR,decreasing=F),]
    #enrich_df<-enrich_df[enrich_df$FDR<0.05,]
    enrich_result_list[[i]]<-enrich_df
}


#将enrich_result_list中的数据合并
enrich_result_df<-do.call(rbind,enrich_result_list)
enrich_result_df$logFDR <- -log10(enrich_result_df$FDR)
enrich_result_df$logFDR[enrich_result_df$logFDR=="Inf"]<-max(enrich_result_df$logFDR[!is.infinite(enrich_result_df$logFDR)])

#输出enrich_result_df
write.csv(enrich_result_df,"D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/enrich_result_df.csv",row.names = F)


outputDirectory <- "D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/gene_enrichment"

enrich_result_list<-list()
for(i in colnames(gene_list)){
  print(i)
    #判断是否有Project_+i文件夹，如果没有，那么跳过
    if(!file.exists(paste0(outputDirectory,"/Project_",i,"/enrichment_results_",i,".txt"))){
      enrich_df<-NULL
    }else{
    #读取Project_+i的文件夹中，enrichment_results_+i+.txt的文件
    enrich_df<-read.delim(paste0(outputDirectory,"/Project_",i,"/enrichment_results_",i,".txt"),header = T)
    enrich_df<-enrich_df[enrich_df$FDR<0.05,]
    enrich_df$cancer<-i
    enrich_df$method<-"pos"
    #
    }
enrich_result_list[[i]]<-enrich_df
}


#将enrich_result_list中的数据合并
enrich_result_df<-do.call(rbind,enrich_result_list)
enrich_result_df$logFDR <- -log10(enrich_result_df$FDR)
enrich_result_df$logFDR[enrich_result_df$logFDR=="Inf"]<-max(enrich_result_df$logFDR[!is.infinite(enrich_result_df$logFDR)])

#输出enrich_result_df
write.csv(enrich_result_df,"D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/enrich_result_df_reactome.csv",row.names = F)