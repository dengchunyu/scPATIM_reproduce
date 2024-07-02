gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')

rolypoly_df<-list()
for(i in gwass){
result_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/rolypoly/{i}_rolypoly.RData"))
load(result_file)
df<-rolypoly_result$bootstrap_results
rownames(df)<-df$annotation
df<-df[c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg"),]
rolypoly_df[[i]]<-df$bp_value

}
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
rolypoly_df<-as.data.frame(rolypoly_df)
colnames(rolypoly_df)<-cancers
rownames(rolypoly_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")

##sclinker
sclinker_df<-list()
for(i in cancers){
    result_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/{i}_cell_types_pvalues.txt"))
    df<-read.table(result_file,header=T)
    #df<-df$Enrichment_p
    df$Cell_Type <- sub(".*_", "", df$Cell_Type)
    rownames(df)<-df$Cell_Type
    print(df$Cell_Type)
    df<-df[c("Plasma","DC" ,"Mono","Mast","NK","T.cyt","Th","Treg"),]
    sclinker_df[[i]]<-df$Enrichment_p
}
sclinker_df<-as.data.frame(sclinker_df)
colnames(sclinker_df)<-cancers
rownames(sclinker_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")

##scPagwas_df
scPagwas_df<-list()

for(i in cancers){
    result_file= file.path(glue::glue('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.2.scPagwas_orignal/{i}/_celltypes_bootstrap_results.csv'))
    df<-read.csv(result_file)
    rownames(df)<-df$annotation
    df<-df[c("B/Plasma","DC" ,"Mac/Mono","Mast","NK","T.cyt","Th","Treg"),]
    scPagwas_df[[i]]<-df$bp_value
}
scPagwas_df<-as.data.frame(scPagwas_df)
colnames(scPagwas_df)<-cancers
rownames(scPagwas_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")


##scDRS
scPatim_df<-list()
scDRS_df<-list()
Fusion_df<-list()

for(i in 1:15){
cancer<-cancers[i]
gwas<-gwass[i]
    result_file= file.path(glue::glue('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/{cancer}/Merged_celltype_pvalue.csv'))
    df<-read.csv(result_file)
    rownames(df)<-df$celltype

result_file= file.path(glue::glue('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/{gwas}_fusion_celltypeP.csv'))
    df3<-read.csv(result_file)
    rownames(df3)<-df$celltype

result_file= file.path(glue::glue('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/scDRS/{gwas}_mergecelltypeP.csv'))
    df2<-read.csv(result_file)
    rownames(df2)<-df$celltype

    df<-df[c("B/Plasma","DC" ,"Mac/Mono","Mast","NK","T.cyt","Th","Treg"),]
    scPatim_df[[i]]<-df$pvalue
    df2<-df2[c("B/Plasma","DC" ,"Mac/Mono","Mast","NK","T.cyt","Th","Treg"),]
    scDRS_df[[i]]<-df2$assoc_mcp

    df3<-df3[c("B/Plasma","DC" ,"Mac/Mono","Mast","NK","T.cyt","Th","Treg"),]
    Fusion_df[[i]]<-df3$assoc_mcp
}
    

scPatim_df<-as.data.frame(scPatim_df)
colnames(scPatim_df)<-cancers
rownames(scPatim_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")

scDRS_df<-as.data.frame(scDRS_df)
colnames(scDRS_df)<-cancers
rownames(scDRS_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")

Fusion_df<-as.data.frame(Fusion_df)
colnames(Fusion_df)<-cancers
rownames(Fusion_df)<-c("B.Plasma","DC" ,"Mac.Mono","Mast","NK","T.cyt","Th","Treg")

####整合结果
scPatim_df$Celltype<-rownames(scPatim_df)
scPatim_df<-reshape::melt(scPatim_df,id.vars="Celltype")
colnames(scPatim_df)<-c("Celltype","Cancers","pvalue")
scPatim_df$Method<-"scPATIM"

scDRS_df$Celltype<-rownames(scDRS_df)
scDRS_df<-reshape::melt(scDRS_df,id.vars="Celltype")
colnames(scDRS_df)<-c("Celltype","Cancers","pvalue")
scDRS_df$Method<-"scDRS"

Fusion_df$Celltype<-rownames(Fusion_df)
Fusion_df<-reshape::melt(Fusion_df,id.vars="Celltype")
colnames(Fusion_df)<-c("Celltype","Cancers","pvalue")
Fusion_df$Method<-"Fusion"

sclinker_df$Celltype<-rownames(sclinker_df)
sclinker_df<-reshape::melt(sclinker_df,id.vars="Celltype")
colnames(sclinker_df)<-c("Celltype","Cancers","pvalue")
sclinker_df$Method<-"sclinker"


rolypoly_df$Celltype<-rownames(rolypoly_df)
rolypoly_df<-reshape::melt(rolypoly_df,id.vars="Celltype")
colnames(rolypoly_df)<-c("Celltype","Cancers","pvalue")
rolypoly_df$Method<-"rolypoly"

scPagwas_df$Celltype<-rownames(scPagwas_df)
scPagwas_df<-reshape::melt(scPagwas_df,id.vars="Celltype")
colnames(scPagwas_df)<-c("Celltype","Cancers","pvalue")
scPagwas_df$Method<-"scPagwas"

mm_df<-rbind(scPagwas_df,rolypoly_df,sclinker_df,Fusion_df,scPatim_df,scDRS_df)
mm_df$logp<- -log10(mm_df$pvalue)
write.csv(mm_df,file="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/cancer_mult_method_presult.csv")
