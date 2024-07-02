library(readr)
library(dplyr)
immsnp<-read_table('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_filter5e2.tsv',col_names = F)

len_snp<-list()
#576870

first_list<-c("finngen_r7_c3_breast_herneg_exallc","finngen_r7_c3_breast_exallc","finngen_r7_c3_bronchus_lung_exallc","finngen_r7_c3_urinary_tract_exallc","finngen_r7_c3_melanoma","finngen_r7_c3_lung_nonsmall_exallc","finngen_r7_c3_colon_exallc","finngen_r7_c3_cervix_uteri_exallc","finngen_r7_c3_colon_adeno_exallc","finngen_r7_c3_bladder_exallc","finngen_r7_cd2_insitu_breast_exallc","finngen_r7_c3_rectum_exallc","finngen_r7_c3_kidney_notrenalpelvis_exallc","finngen_r7_c3_thyroid_gland_exallc","finngen_r7_c3_ovary_exallc","finngen_r7_c3_nsclc_adeno_exallc","finngen_r7_c3_pancreas_exallc","finngen_r7_c3_stomach_exallc","finngen_r7_c3_nsclc_squam_exallc","finngen_r7_c3_anus_analcanal_exallc","finngen_r7_c3_baseoftongue_exallc","finngen_r7_c3_biliary_tract_exallc","finngen_r7_c3_brain_exallc","finngen_r7_c3_corpus_uteri_exallc","finngen_r7_c3_digestive_organs_exallc","finngen_r7_c3_endocrine_exallc","finngen_r7_c3_eye_adnexa_exallc","finngen_r7_c3_eye_brain_neuro_exallc","finngen_r7_c3_floorofmouth_exallc","finngen_r7_c3_gbm_exallc","finngen_r7_c3_heart_mediastinum_pleura_exallc","finngen_r7_c3_larynx_exallc","finngen_r7_c3_lip_exallc","finngen_r7_c3_lip_oral_pharynx_exallc","finngen_r7_c3_liver_intrahepatic_bile_ducts_exallc","finngen_r7_c3_male_genital_exallc","finngen_r7_c3_meninges_exallc","finngen_r7_c3_mestothel_softtissue_exallc","finngen_r7_c3_mouthnas_exallc","finngen_r7_c3_oesophagus_exallc","finngen_r7_c3_parotidgland_exallc","finngen_r7_c3_prostate_exallc","finngen_r7_c3_renal_pelvis_exallc","finngen_r7_c3_skin_exallc","finngen_r7_c3_small_intestine_exallc","finngen_r7_c3_tongue_exallc","finngen_r7_c3_tonguenas_exallc","finngen_r7_c3_tonsil_exallc","finngen_r7_c3_vagina_exallc","finngen_r7_cd2_benign_brain_cns_exallc","finngen_r7_cd2_benign_brain_supratent_exallc","finngen_r7_cd2_benign_eye_exallc","finngen_r7_cd2_eye_nos_exallc")


for (i in first_list){
	gwas<- read_table(paste0("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/",i,"_prune.txt"))
	snps<-intersect(gwas$rsid,immsnp$X2)
	len_snp[i]<-length(snps)
    print(length(snps))
    if(length(snps)>100000){
        gwas<-gwas[gwas$rsid %in% snps,]
    write.table(gwas,file= paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/",i,"_prune_isnp.txt"),row.names=F,quote=F)	
    }
}


second_list<-c('bbj-a-158','bbj-a-117','ieu-b-4810','ieu-b-4965','bbj-a-107','ieu-b-4969','ieu-b-4954','bbj-a-133','ieu-a-1013','ieu-b-4963','GCST007728','ieu-b-4874','bbj-a-119','ieu-a-822','ukb-b-1316','ebi-a-GCST006464','ieu-b-4912','bbj-a-140')
len_snp<-list()
for (i in second_list){

	gwas<- read_table(paste0("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/",i,"_prune.txt"))
    #'GCST007728'
    #gwas<- read_table(paste0("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/",i,".txt"))
	snps<-intersect(gwas$rsid,immsnp$X2)
	#len_snp[i]<-length(snps)
    print(length(snps))
    if(length(snps)>100000){
        gwas<-gwas[gwas$rsid %in% snps,]
    write.table(gwas,file= paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/",i,"_prune_isnp.txt"),row.names=F,quote=F)
    }
}