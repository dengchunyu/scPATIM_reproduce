# 画泛癌不同MP的影响的热图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(stringr)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(stringr)
cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass<-c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')

# 添加数据到数据框
Full_Name <- c('DC', 'B cells', 'Regulatory T cells', 'Germinal center B cells', 'Tem/Temra cytotoxic T cells', 'CD16- NK cells', 'Proliferative germinal center B cells', 'Tcm/Naive cytotoxic T cells', 'DC2', 'pDC', 'Alveolar macrophages', 'Non-classical monocytes', 'Tcm/Naive helper T cells', 'Tem/Trm cytotoxic T cells', 'NK cells', 'Trm cytotoxic T cells', 'Type 17 helper T cells', 'Tem/Effector helper T cells', 'CD16+ NK cells', 'Migratory DCs', 'Naive B cells', 'Mast cells', 'Macrophages', 'Classical monocytes', 'ILC', 'Plasma cells', 'Memory B cells', 'Intermediate macrophages', 'gamma-delta T cells', 'MAIT cells', 'Erythrophagocytic macrophages','HSC/MPP')
Abbreviation <- c('DC', 'Bcells' , 'Treg', 'Bgc', 'Tem/Temra.cyt', 'NK.CD16neg', 'Bpgc', 'Tcm/N.cyt', 'DC2', 'pDC', 'Mac', 'Mono.nc', 'Tcm/N.h', 'Tem/Trm.cyt', 'NK', 'Trm.cyt', 'Th17', 'Tem/Eff.h', 'NK.CD16pos', 'DC.mig', 'Bnaive', 'Mast', 'Mac', 'Mono.c', 'ILC', 'Plasma', 'Bmem', 'Mac.inter', 'Tgd','MAIT','Ery.Mac','HSC/MPP')
merge_Abbreviation <- c('DC', 'B/Plasma' , 'Treg', 'B/Plasma', 'T.cyt', 'NK', 'B/Plasma', 'T.cyt', 'DC', 'DC', 'Mac/Mono', 'Mac/Mono', 'Th', 'T.cyt', 'NK', 'T.cyt', 'Th', 'Th', 'NK', 'DC', 'B/Plasma', 'Mast', 'Mac/Mono', 'Mac/Mono', 'ILC', 'B/Plasma', 'B/Plasma', 'Mac/Mono', 'Tgd','MAIT','Mac/Mono','HSC/MPP')
cl_annotation <- data.frame(
  Abbreviation = Abbreviation,
  Full_Name = Full_Name,
  merge_Abbreviation=merge_Abbreviation,
  stringsAsFactors = FALSE
)
zscoref<-function(vec){
mean_val <- mean(vec)
sd_val <- sd(vec)
zscore_vec <- (vec - mean_val) / sd_val   
}

setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis")
mp_sumscore<-list()
anno_list<-list()
for(i in 1:14){
    cancer<-cancers[i]
    gwas<-gwass[i]
    print(cancer)
    load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_pagwas.Rdata"))

    Pagwas$Celltype_anno$scPagwas.gPAS.score<-Pagwas$scPagwas.gPAS.score
    Celltype_anno<-Pagwas$Celltype_anno

    df<-read.csv(paste0(cancer,"_merge_average_pscores.csv"))
    if(all(unique(Celltype_anno$annotation) %in% cl_annotation$Full_Name)==TRUE){
        for(j in unique(Celltype_anno$annotation)){
        print(j)
        Celltype_anno$annotation[Celltype_anno$annotation==j]<-cl_annotation$Abbreviation[cl_annotation$Full_Name==j]
        Celltype_anno$merge_anno<-Celltype_anno$annotation
        }
        for(j in unique(Celltype_anno$annotation)){
        Celltype_anno$merge_anno[Celltype_anno$merge_anno==j]<-cl_annotation$merge_Abbreviation[cl_annotation$Abbreviation==j]
        }
    }else{
        stop("some cells no include!")
    }
    save(Celltype_anno,file=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_Celltype_anno.Rdata"))

    sig<- df$celltype[which(df$pvalue<0.05)]
    sig_annocell<-Celltype_anno[Celltype_anno$merge_anno %in% sig,]
    mp_sigcell<-Pagwas$Pathway_single_results[,sig_annocell$cellnames]
    
    top_sig_annocell<-sig_annocell[order(sig_annocell$scPagwas.gPAS.score,decreasing=T)[1:500],]
    top_sig_annocell$cancer<-cancer
    mat_top<-Pagwas$Pathway_single_results[,top_sig_annocell$cellnames]

    mat_top_processed<-apply(mat_top,2, zscoref)
    mp_sumscore[[i]]<-mat_top_processed
    anno_list[[i]]<-top_sig_annocell
}

mp_list<-lapply(mp_sumscore,function(x) rownames(x))
names(mp_score_list)<-cancers
names(anno_list)<-cancers

mp_list<-Reduce(intersect,mp_list)
mp_score_list<-lapply(mp_sumscore,function(x) x[mp_list,])

mp_mat<-Reduce(cbind,mp_score_list)
anno_mat<-Reduce(rbind,anno_list)
save(mp_score_list,anno_list,file="mp_anno_list.RData")
save(mp_mat,anno_mat,file="mp_anno_topmat.RData")

######画图
setwd("D:/Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
load("D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_anno_topmat.RData")
heatmap_legend_param <- list(
  title = "Hereditary effect",
  at = seq(min(mp_mat), max(mp_mat), length.out = 5),  # 设置刻度
  limits = c(min(mp_mat), max(mp_mat)),
  labels =round(seq(min(mp_mat), max(mp_mat), length.out = 5),3)  # 设置数值范围
)

ha = HeatmapAnnotation(
    Celltype =anno_mat$annotation ,
    Sub_celltype=anno_mat$merge_anno,
    cancers=anno_mat$cancer
    )

pdf("heatmap_merge_pathway.pdf",width=16,height=13)
print(ComplexHeatmap::Heatmap(as.matrix(mp_mat),name ="Hereditary effect",col=c("#116979", "#DEE3E2", "#DE7119"),heatmap_legend_param = heatmap_legend_param,border = NA,
        column_title=NULL,show_column_names=FALSE,cluster_columns=TRUE,
        top_annotation = ha,
        show_column_dend=TRUE))
dev.off()