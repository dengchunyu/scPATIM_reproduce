library(plyr)
library(stringr)
library(Seurat)
library(Matrix)
library(scPagwas)
cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','LiverCancer','ColorectalCancer','EsophagealCancer','GastricCancer','EndometrialCarcinoma')
gwass<- c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','finngen_r7_c3_stomach_exallc','ebi-a-GCST006464')
index<-data.frame(cancers,gwass)
#setwd("D:/OneDrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis")

for(i in 1:15){
    cancer<-cancers[i]
    gwas<-gwass[i]
    print(cancer)
    load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_pagwas.Rdata"))
    Single_data<-CreateSeuratObject(
      Pagwas$data_mat,
      project = "CreateSeuratObject",
      assay = "RNA",
      names.field = 1,
      names.delim = "_")
    scPagwas_pathway <- SeuratObject::CreateAssayObject(data = Pagwas$Pathway_single_results)
    Single_data[["scPagwasPaHeritability"]] <- scPagwas_pathway
    gpas<-Pagwas$scPagwas.gPAS.score
    n<- floor(0.05 * length(gpas))
    
    top_sig_cell<-names(gpas)[order(gpas,decreasing=T)[1:n]]
    Single_data$sig_cell<-"NonSig"
    Single_data$sig_cell[colnames(Single_data) %in% top_sig_cell]<-"Sig"

    Idents(Single_data)<-Single_data$sig_cell
    sig_marker<-FindAllMarkers(Single_data,
      assay = "scPagwasPaHeritability",
      logfc.threshold = 0,
      min.pct=0,
      return.thresh = 1,
      test.use = "wilcox",
      slot = "data")
    sig_marker2<-sig_marker[sig_marker$cluster=="Sig",]
    sig_marker2$cancers<-cancer
    sig_marker2$adjp<-p.adjust(sig_marker2$p_val,method="BH")
save(sig_marker2,file=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/",cancer,"_",gwas,"_sigcell_marker_PaHeri.RData"))
}

#整合差异分析的结果，输出top5的mp
cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','LiverCancer','ColorectalCancer','EsophagealCancer','GastricCancer','EndometrialCarcinoma')
library(dplyr)
mp_list<-list()
for(i in 1:15){
    cancer<-cancers[i]
    gwas<-gwass[i]
    print(cancer)
    load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/",cancer,"_",gwas,"_sigcell_marker_PaHeri.RData"))
    top5_genes <- sig_marker2 %>%
  top_n(5, wt =avg_log2FC) %>%
  ungroup()
    mp_list[[i]]<-top5_genes
}
mp_df<- Reduce(rbind,mp_list)
write.csv(mp_df,file="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/Cancers_top_HCPs.csv")