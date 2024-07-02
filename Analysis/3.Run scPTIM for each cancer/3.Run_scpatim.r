library(Seurat)
library(Matrix)
library(scPagwas)
Args <- commandArgs(T)
num = print(Args[1])
#'HeadandNeck',
cancers=c('EndometrialCarcinoma','KidneyCancer','ProstateCancer','Pancreatic','BreastCancer','EsophagealCancer','ThyroidCancer','BladderCancer','ColorectalCancer','LiverCancer','LungCancer','Melanoma','OvarianCancer','GastricCancer')
#cancers=c('BladderCancer','ColorectalCancer','GastricCancer')
#gwass=c('finngen_r7_c3_bladder_exallc','finngen_r7_c3_colorectal','finngen_r7_c3_stomach_exallc')
gwass=c('ebi-a-GCST006464','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_breast_exallc','ieu-b-4912','bbj-a-117','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_colorectal','finngen_r7_c3_liver_intrahepatic_bile_ducts_exallc','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_ovary_exallc','finngen_r7_c3_stomach_exallc')
##第二轮计算

gwass=c("ieu-b-4810","ieu-b-4969","ukb-d-C3_MELANOMA_SKIN","ieu-b-4954","bbj-a-133","ieu-b-4874","ukb-b-1316","ieu-b-4963","bbj-a-119","ieu-a-822","bbj-a-158","ieu-b-4965","bbj-a-107","finngen_r7_c3_oesophagus_exallc",)#"ukb-b-18843",
cancers=c("BreastCancer","Melanoma","Melanoma","LungCancer","LungCancer","BladderCancer","KidneyCancer","OvarianCancer","GastricCancer","GastricCancer","LiverCancer","ColorectalCancer","ColorectalCancer","EsophagealCancer")#"ProstateCancer"

gwass=c("GCST90041886","bbj-a-102","bbj-a-148","GCST90018929","bbj-a-139","bbj-a-140","GCST90018838","ukb-b-18843")
cancers=c("BreastCancer","BreastCancer","ProstateCancer","ThyroidCancer","OvarianCancer","Pancreatic","EndometrialCarcinoma","ProstateCancer")

#setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.scPagwas_run/")
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.2.scPagwas_orignal/")
num<-as.numeric(num)
gwas<-gwass[num]
cancer<-cancers[num]
gwas_path <- file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/1.Gwas_input/{gwas}_gwas_pagwas.RData"))
scmat_path <- file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/{cancer}_single.pagwas.RData"))

load(gwas_path)
Pagwas1<-Pagwas
rm(Pagwas)

load(scmat_path)
Pagwas<- c(Pagwas1,Pagwas)

rm(Pagwas1)
gc()


output.prefix=gwas
output.dirs=cancer
if (!dir.exists(output.dirs)) {
    dir.create(output.dirs)
}
if (!dir.exists(output.dirs)) {
  dir.create(paste0("./", output.dirs, "/temp"))
}

if(gwas %in% c('ebi-a-GCST006464','ieu-b-4912','bbj-a-117','ieu-b-4810','ieu-b-4969')){
    Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation_hg37
    )
    }else{
        #hg38
Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation
    )
    }


  
Pagwas <- Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      singlecell = T,
      celltype = F,
      backingpath=paste0("./", output.dirs, "/temp"),
      n.cores=1
      )

    #Pagwas$Pathway_ld_gwas_data <- NULL
    Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )
save(Pagwas,file=paste0("./", output.dirs, "/",output.prefix, "_pagwas.Rdata"))

  utils::write.table(Pagwas$Pathway_sclm_results,
    file = paste0(
      "./", output.dirs, "/",
      output.prefix, "_Pathway_singlecell_lm_results.txt"
    ),
    quote = F, sep = "\t"
  )