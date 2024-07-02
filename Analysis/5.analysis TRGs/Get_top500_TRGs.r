setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')

gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
gene_list<-list()
for(i in 1:14){
    cancer<-cancers[i]
    gwas<-gwass[i]
  PCC<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_gene_PCC.csv"))
  gene_list[[i]]<-PCC$X[order(PCC$weight_pcc, decreasing =T)[1:500]]
}

gene_list<-as.data.frame(gene_list)
colnames(gene_list)<-cancers
write.csv(gene_list,file="cancers_top500.csv")
save(gene_list,file="cancers_top500.RData")