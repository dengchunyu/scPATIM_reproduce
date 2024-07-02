## LDSC+sclinker计算不同细胞类型的pvalue

获得基因区域数据

```R
library(tidyverse)
library(data.table)
library(dplyr)
##10kb
gene_coordinates0 <-
  read_tsv("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-50000<0,0,X3-50000),end=X4+50000) %>%
  select(X2,start,end,6,1) %>%
  rename(chr="X2", Gene="X6",EntrezGene="X1") %>%
  mutate(chr=paste0("chr",chr))

##abc model
df_pre = data.frame(fread("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz"))
df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
unique(df_pre$CellType)
tissuename ="BLD"
tissuenames2 = as.character(read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/ABC.listbloodQC.txt", header=F)[,1])
tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
df = df_pre[tissue_ids, ]
df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
colnames(df2) = c("chr", "start", "end", "Gene")

##enhancer model
Enhancer = read.table("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/Roadmap_Enhancers_Blood.txt",header=F)
df3<-Enhancer[,1:4]
colnames(df3)<-c("chr","start","end","Gene")
df3 = rbind(df3,df2)
gene_coordinates<-df3

#将gene_coordinates0中的EntrezGene列加入到gene_coordinates中
gene_coordinates<-merge(gene_coordinates,gene_coordinates0[,c("Gene","EntrezGene")],by="Gene",all.x=TRUE)
#删除na
gene_coordinates<-gene_coordinates[!is.na(gene_coordinates$EntrezGene),]
save(gene_coordinates,file="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/ABC_gene_coordinates.RData")
#输出gene_coordinates到csv
write.table(gene_coordinates,file="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/ABC_gene_coordinates.csv",sep="\t",quote=F,row.names=F,col.names=T)
```

计算每个基因的平均表达量，以及输出bed

```python
import pandas as pd
import numpy as np
import os
import scanpy as sc
import gc
import numpy as np
import anndata
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#
cancers=['HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']
#,
gene_coordinates=pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/ABC_gene_coordinates.csv",sep="\t")

for i in cancers:
    print(i)
    results_file=os.path.join(BaseDirectory,i+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    unique_categories = np.unique(adata.obs['merge_celltype_annotation'])
    average_expression_dict = {}
    for category in unique_categories:
        cells_in_category = adata.obs['merge_celltype_annotation'] == category
        average_expression = np.mean(adata.X[cells_in_category], axis=0)
        if np.ndim(average_expression) > 1:
            average_expression = average_expression.A1
        average_expression_dict[category] = average_expression
    average_expression_df = pd.DataFrame.from_dict(average_expression_dict, orient='columns')
    average_expression_df.index=adata.var_names
    bedfile="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/"+i
    if not os.path.exists(bedfile):
        os.makedirs(bedfile)
    ave_func(average_expression_df, bedfile)


def ave_func(average_expression, path):
    average_expression = average_expression.dropna()
    average_expression = average_expression.loc[average_expression.sum(axis=1) != 0]
    top10_function(average_expression, path)

def top10_function(exp, path):
    exp['Gene'] = exp.index
    exp = exp.groupby('Gene').filter(lambda x: len(x) == 1)
    exp = exp.melt(id_vars='Gene', var_name='column', value_name='Expr')
    exp = exp.assign(Expr_sum_mean=exp['Expr'] * 1e6 / exp.groupby('column')['Expr'].transform('sum'))
    exp = exp.assign(specificity=exp['Expr_sum_mean'] / exp.groupby('Gene')['Expr_sum_mean'].transform('sum'))
    exp2 = exp.merge(gene_coordinates, on='Gene', how='inner')
    n_genes = exp2['EntrezGene'].nunique()
    n_genes_to_keep = round(n_genes * 0.1)
    filtered_exp = exp2[exp2['Expr_sum_mean'] > 1]
    ldsc_bedfile(filtered_exp, 'column', n_genes_to_keep, path)
    print("success!")

def write_group(df, Cell_type, path):
    df = df[['column', 'chr', 'start', 'end', 'EntrezGene']]
    os.makedirs(path, exist_ok=True)
    file_path = os.path.join(path, f"{make_valid_filename(Cell_type)}.bed")
    df.iloc[:, 1:].to_csv(file_path, sep='\t', header=False, index=False)
    return df

def ldsc_bedfile(d, Cell_type, n_genes_to_keep, path):
    d_spe = d.groupby(Cell_type).apply(lambda x: x.nlargest(n_genes_to_keep, 'specificity')).reset_index(drop=True)
    d_spe.groupby(Cell_type).apply(lambda x: write_group(x, x.name, path))

def make_valid_filename(filename):
    keepcharacters = (' ', '.', '_')
    return "".join(c if c.isalnum() or c in keepcharacters else "_" for c in filename)


```

### 3.5 get annot files

```shell
cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/
source activate mypy
bunzip2 w_hm3.snplist.bz2
tail -n +2 w_hm3.snplist | cut -f 1 > hm_snp.txt

cd /share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/1000G_EUR_Phase3_baseline/
#gunzip baseline.*.annot.gz
for file in baseline.*.annot
do
awk 'NR > 1{print "chr"$1"\t"$2"\t"$2"\t"$3}' $file >> tmp.bed
done
sortBed -i tmp.bed > 1000genomes_phase3_SNPs.bed2
rm tmp.bed

cancers=['HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']


#!/bin/bash
#SBATCH -e EsophagealCancer.err
#SBATCH -o EsophagealCancer.out
#SBATCH -J EsophagealCancer
#SBATCH -w in002
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
#Mac_Mono Mast Treg Th
cancer=EsophagealCancer
#f=T.cyt
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/${cancer}
path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
all_snps="1000G_EUR_Phase3_baseline/1000genomes_phase3_SNPs.bed2"
all_annotations="1000G_EUR_Phase3_baseline/"
plink_file="1000G_EUR_Phase3_plink/"
hapmap_snps="hm_snp.txt"
weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
frq="1000G_Phase3_frq/1000G.EUR.QC."
echo $cancer

for f in *.bed
do
echo $f
intersectBed -c -a $path_name$all_snps -b $f > $f".1000genomes.intersect"
awk '{if($5!=0) print $4}' $f".1000genomes.intersect" > $f".1000genomes.intersect.snp"
mkdir $f"_tissue_dir"
rm $f".1000genomes.intersect"
cd $f"_tissue_dir"
for j in $path_name$all_annotations/*.annot
do
echo $j
file_name=`basename $j`
perl $path_name/fast_match2_minimal.pl ../$f".1000genomes.intersect.snp" $f $j > $file_name
done
gzip *annot
for i in {1..22}
do
/share/pub/dengcy/software/ldsc-master/ldsc.py --l2 --bfile $path_name$plink_file/1000G.EUR.QC.$i --ld-wind-cm 1 --print-snps $path_name$hapmap_snps --annot baseline.$i.annot.gz --out ./baseline.$i
done
cd ..
rm $f".1000genomes.intersect.snp"
done


```

### 3.6 ldsc

```shell
#!/bin/bash
#SBATCH -e Tcells.err
#SBATCH -o Tcells.out
#SBATCH -J Tcells
#SBATCH -w in010
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate ldsc
#'HeadandNeck' 
cancers=('BreastCancer' 'ProstateCancer' 'Melanoma' 'LungCancer')
cancers=('BladderCancer' 'KidneyCancer' 'ThyroidCancer' 'OvarianCancer')
cancers=('Pancreatic' 'GastricCancer' 'LiverCancer' 'ColorectalCancer' 'EsophagealCancer' 'EndometrialCarcinoma')
cancers=('EndometrialCarcinoma')
#'ieu-b-4912' 
gwass=('finngen_r7_c3_breast_exallc' 'finngen_r7_c3_prostate_exallccd ' 'finngen_r7_c3_melanoma' 'finngen_r7_c3_lung_nonsmall_exallc')
gwass=('finngen_r7_c3_bladder_exallc' 'finngen_r7_c3_kidney_notrenalpelvis_exallc' 'finngen_r7_c3_thyroid_gland_exallc' 'finngen_r7_c3_ovary_exallc')
gwass=( 'finngen_r7_c3_pancreas_exallc' 'finngen_r7_c3_stomach_exallc' 'bbj-a-158' 'finngen_r7_c3_colorectal' 'bbj-a-117')
gwass=('ebi-a-GCST006464')
#nodes=("in002" "in003" "in004")
#node_count=${#nodes[@]}
#counter=0
for ((i=0; i<${#cancers[@]}; i++))
do
  cancer="${cancers[i]}"
  gwas="${gwass[i]}"

  sumstats="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/heritability/${gwas}_isnp.sumstats.gz"
  path_name="/share/pub/dengcy/GWAS_Multiomics/Revise_2023.4.11/sclinker/processed_data/"
  weights="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
  frq="1000G_Phase3_frq/1000G.EUR.QC."
  all_annotations="1000G_EUR_Phase3_baseline"
  cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/${cancer}

  for f in *_tissue_dir
  do
  echo $f
  gwas_name=`basename $sumstats | cut -d "." -f 1`
  echo $gwas_name
  cd $f
  /share/pub/dengcy/software/ldsc-master/ldsc.py --h2 $sumstats --ref-ld-chr $path_name$all_annotations/baseline.,baseline. --w-ld-chr $path_name$weights --overlap-annot --frqfile-chr $path_name$frq --print-coefficients --out ./${gwas_name}_${f}
  cd ..
  done
done
  mkdir log_get_pvalues

```

```R
library("tidyverse")
library("stringr")
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
for(cancer in cancers){
 setwd(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/",cancer))
files <- list.files(".",pattern=".tissue_dir.results",full.names = TRUE,recursive=T)
d <- data_frame(filename=files) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames=gsub(".bed_tissue_dir.results","",basename(filename)),
           makenames=gsub(".bed_continuous_tissue_dir.results","",basename(makenames))) %>% unnest()
  d <- d %>%  filter(Category=="L2_1") %>% mutate(P=1-pnorm(`Coefficient_z-score`)) %>%
  mutate(Trait=sub("_.*","",makenames),Cell_Type=gsub("^_","",str_extract(makenames, "_.*"))) %>%
  select(Trait,Cell_Type,Enrichment,Enrichment_std_error,Enrichment_p,P) %>% arrange(Enrichment_p)
  d$cancer<-cancer
  d$Cell_Type<-sub("isnp_", "",  d$Cell_Type)
write_tsv(d,path=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/sclinker/",cancer,"_cell_types_pvalues.txt"))
 result_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/rolypoly/{gwas}_rolypoly.RData"))

}
```
