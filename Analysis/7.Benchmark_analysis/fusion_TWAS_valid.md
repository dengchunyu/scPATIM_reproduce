## TWAS analysis

流程学习地址：
http://gusevlab.org/projects/fusion/#installation

### create sumstats file
```r
install.packages('/share/pub/dengcy/software/fusion_twas-master/plink2R-master/plink2R/',repos=NULL)

rm(list=ls())
library(data.table)
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')

gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')

for(gwas in gwass){
	gwasfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/{gwas}.txt"))
	gwasfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/{gwas}.txt"))

gwas_data <- bigreadr::fread2(gwasfile)
gwas_data$beta<-as.numeric(gwas_data$beta)
gwas_data$se<-as.numeric(gwas_data$se)
gwas_data$maf<-as.numeric(gwas_data$maf)
gwas_data$pos<-as.numeric(gwas_data$pos)
gwas_data$Z<-gwas_data$beta/gwas_data$se
gwas_data1<-gwas_data[,c("rsid","REF", "ALT","Z")]
colnames(gwas_data1)<-c("SNP","A1","A2","Z")
resultfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/{gwas}.sumstats"))
write.table(gwas_data1,file=resultfile,sep=" ",quote=F)
}

``````
### run TWAS

```shell
#PBS -N fusion
#PBS -q workq
#PBS -l mem=150gb
#PBS -l ncpus=2


nodes=("in009" "in004" "in003")
gwass=('ieu-b-4912' 'finngen_r7_c3_breast_exallc' 'finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma' 'finngen_r7_c3_lung_nonsmall_exallc' 'finngen_r7_c3_bladder_exallc' 'finngen_r7_c3_kidney_notrenalpelvis_exallc' 'finngen_r7_c3_thyroid_gland_exallc' 'finngen_r7_c3_ovary_exallc' 'finngen_r7_c3_pancreas_exallc' 'finngen_r7_c3_stomach_exallc' 'bbj-a-158' 'finngen_r7_c3_colorectal' 'bbj-a-117' 'ebi-a-GCST006464')
#!/bin/bash
#SBATCH -e twas14.err
#SBATCH -o twas14.out
#SBATCH -J twas14
#SBATCH -w in003
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion
source activate myr
gwas=ebi-a-GCST006464
for chr in $(seq 1 22);
do
Rscript /share/pub/dengcy/software/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/${gwas}.sumstats \
--weights  ./WEIGHTS/GTEx.Whole_Blood.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr /share/pub/dengcy/software/fusion_twas-master/LDREF/1000G.EUR. \
--chr $chr \
--out /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/${gwas}.$chr.dat;
done
```

### 整合结果

```R

gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
for(gwas in gwass){
	re_list<-list()
	for(i in 1:22){
		re<-read.table(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/",gwas,".",i,".dat"),header=T)
		re_list[[i]]<-re[,c("ID","TWAS.P")]
	}
	re_df<-Reduce(rbind,re_list)
	re_df<-re_df[!duplicated(re_df$ID),]
	re_df$logp<- -log10(re_df$TWAS.P)
	re_df<-re_df[,c("ID","logp")]
	colnames(re_df)<-c("gene_symbol","logp")
	write.csv(re_df,file=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/",gwas,"_fusion_gene.csv"),quote=F)
}

```
计算细胞类型pvalue

```bash
#!/bin/bash
#SBATCH -e twas.err
#SBATCH -o twas.out
#SBATCH -J twas
#SBATCH -w in003
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00

source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion

# 定义数组
cancers=('HeadandNeck' 'BreastCancer' 'ProstateCancer' 'Melanoma' 'LungCancer' 'BladderCancer' 'KidneyCancer' 'ThyroidCancer' 'OvarianCancer' 'Pancreatic' 'GastricCancer' 'LiverCancer' 'ColorectalCancer' 'EsophagealCancer' 'EndometrialCarcinoma')
gwass=('ieu-b-4912' 'finngen_r7_c3_breast_exallc' 'finngen_r7_c3_prostate_exallc' 'finngen_r7_c3_melanoma' 'finngen_r7_c3_lung_nonsmall_exallc' 'finngen_r7_c3_bladder_exallc' 'finngen_r7_c3_kidney_notrenalpelvis_exallc' 'finngen_r7_c3_thyroid_gland_exallc' 'finngen_r7_c3_ovary_exallc' 'finngen_r7_c3_pancreas_exallc' 'finngen_r7_c3_stomach_exallc' 'bbj-a-158' 'finngen_r7_c3_colorectal' 'bbj-a-117' 'ebi-a-GCST006464')

nodes=("in005" "in006" "in004")
node_count=${#nodes[@]}
counter=0

# 循环迭代每个癌症
for ((i=0; i<${#cancers[@]}; i++))
do
  cancer="${cancers[i]}"
  gwas="${gwass[i]}"
  node="${nodes[$counter % $node_count]}"
  # 创建任务脚本文件
  #script_file="${cancer}_task.sh"
	filename="script_${j}.sh"
	#cancer=EndometrialCarcinoma
  #gwas=ebi-a-GCST006464
    # 创建.sh文件并写入内容
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e twas_${cancer}.err
#SBATCH -o twas_${cancer}.out
#SBATCH -J twas_${cancer}
#SBATCH -w ${node}
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
python /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/${cancer}_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/${gwas}_fusion_gene.csv --top_gene_num 1000 --n_ctrl 200 --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/${gwas}_fusion_score.csv --group merge_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/${gwas}_fusion_celltypeP.csv
EOF
  # 提交任务脚本
  sbatch "$filename"
  counter=$((counter + 1))
done

```

### 计算细胞类型显著性

```python
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import pyarrow as pa
import os
import scipy.sparse
import csv
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import norm

BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
cancers=['HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']

gwass=['ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464']
for i in range(0,15):
    print(i)
    gwas=gwass[i]
    cancer=cancers[i]
    results_file=os.path.join(BaseDirectory,cancer+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    drs_df=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/" +gwas+ "_fusion_score.csv")
    drs_df=drs_df.iloc[:,[1,3]]
    drs_df.columns='scDRS_'+drs_df.columns
    drs_df.index=adata.obs.index
    drs_df["merge_celltype_annotation"]=adata.obs["merge_celltype_annotation"]
    ct_value=merge_celltype_p(drs_df['scDRS_pval'], drs_df['merge_celltype_annotation'])
    ct_value.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/Fusion/result/'+cancer+'_Merged_pvalue.csv')

def p_merge(pvalues):
    zvalues = -np.sqrt(2) * norm.ppf(pvalues / 2)
    ztotal = np.mean(zvalues)
    p_total = norm.cdf(-abs(ztotal))
    return p_total

def merge_celltype_p(single_p, celltype):
    celltype_p = pd.DataFrame({'celltype': celltype, 'pvalue': single_p})
    celltype_p = celltype_p.groupby('celltype')['pvalue'].agg(p_merge).reset_index()
    return celltype_p