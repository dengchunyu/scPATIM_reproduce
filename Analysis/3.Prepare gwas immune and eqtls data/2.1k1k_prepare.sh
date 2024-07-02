#############
#1.read the 1k1k data,filter the pvalue
################
#shell

sed 's/ //g' onek1k_eqtl_dataset.tsv > /share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_dataset.tsv

awk '{print $2,$3,$5,$7,$8,$15}' onek1k_eqtl_dataset.tsv > /share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_filterpvalue.tsv
#64374851
cd /share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/
awk '{if($6<=5e-2) print $1,$2,$3,$4,$5}' onek1k_eqtl_filterpvalue.tsv >onek1k_eqtl_filter5e2.tsv

wc -l onek1k_eqtl_filter5e2.tsv
#4077555