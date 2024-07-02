########################################
#finngen数据的数据计算
#######################################
library(scPagwas)
#devtools::install_github("sulab-wmu/scPagwas")
#install.packages("/share/pub/dengcy/software/scPagwas_1.10.3.tar.gz",repos=NULL,type="source")
library(Seurat)
#Args <- commandArgs(T)
#gwas = print(Args[1])
#gwas = "finngen_r7_c3_breast_exallc"
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/1.Gwas_input/")
#output.dirs<-gwas
first_list<-c('finngen_r7_c3_breast_her2neg_exallc' ,'finngen_r7_c3_male_genital_exallc' ,'finngen_r7_c3_prostate_exallc' ,'finngen_r7_c3_breast_herplus_exallc','finngen_r7_c3_urinary_tract_exallc' ,'finngen_r7_c3_melanoma' ,'finngen_r7_c3_lung_nonsmall_exallc' ,'finngen_r7_c3_colon_exallc' ,'finngen_r7_c3_cervix_uteri_exallc' ,'finngen_r7_c3_bladder_exallc' ,'finngen_r7_c3_rectum_exallc' ,'finngen_r7_c3_kidney_notrenalpelvis_exallc' ,'finngen_r7_c3_endocrine_exallc' ,'finngen_r7_c3_thyroid_gland_exallc' ,'finngen_r7_c3_lip_oral_pharynx_exallc' ,'finngen_r7_c3_ovary_exallc' ,'finngen_r7_c3_pancreas_exallc' ,'finngen_r7_c3_stomach_exallc' ,'finngen_r7_c3_meninges_exallc' ,'finngen_r7_c3_mestothel_softtissue_exallc' ,'finngen_r7_c3_brain_exallc' ,'finngen_r7_c3_liver_intrahepatic_bile_ducts_exallc' ,'finngen_r7_c3_small_intestine_exallc' ,'finngen_r7_c3_oesophagus_exallc' ,'finngen_r7_c3_eye_adnexa_exallc','finngen_r7_c3_breast_herneg_exallc',
'finngen_r7_c3_breast_exallc',
'finngen_r7_c3_colorectal',
'finngen_r7_c3_nsclc_squam_exallc',
'finngen_r7_c3_nsclc_adeno_exallc')



lapply(first_list,function(gwas){
gwasfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/{gwas}_isnp.txt"))
out_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/1.Gwas_input/{gwas}_gwas_pagwas.RData"))
Pagwas<-list()
gwas_data <- bigreadr::fread2(gwasfile)
colnames(gwas_data)<-gwas_data[2,]
gwas_data<-gwas_data[-c(1,2),]
gwas_data$beta<-as.numeric(gwas_data$beta)
gwas_data$se<-as.numeric(gwas_data$se)
gwas_data$maf<-as.numeric(gwas_data$maf)
gwas_data$pos<-as.numeric(gwas_data$pos)
write.table(gwas_data,file=gwasfile)

Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data
  )
Pagwas$snp_gene_df <- SnpToGene(gwas_data = Pagwas$gwas_data,
        block_annotation = block_annotation,
        marg = 10000)
save(Pagwas,file=out_file)
   })

second_list<-c('bbj-a-158','bbj-a-117','ieu-b-4810','ieu-b-4965','bbj-a-107','ieu-b-4969','ieu-b-4954','bbj-a-133','ieu-a-1013','ieu-b-4963','ieu-b-4874','bbj-a-119','ieu-a-822','ukb-b-1316','ebi-a-GCST006464','ieu-b-4912','bbj-a-140')

lapply(second_list,function(gwas){
gwasfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/{gwas}_prune_isnp.txt"))
out_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/1.Gwas_input/{gwas}_gwas_pagwas.RData"))
Pagwas<-list()
gwas_data <- bigreadr::fread2(gwasfile)
gwas_data$maf<-0.1
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data
  )
Pagwas$snp_gene_df <- SnpToGene(gwas_data = Pagwas$gwas_data,
        block_annotation = block_annotation,
        marg = 10000)
save(Pagwas,file=out_file)
   })
#install.packages(, repos = NULL, type = 'source')

