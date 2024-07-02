## rolypoly valid

### 前期数据准备
计算平均表达数据

计算每个细胞类型的平均表达

```python
import numpy as np
import anndata
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#'HeadandNeck','BreastCancer','ProstateCancer',
cancers=['Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']
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
    avg_file=os.path.join(BaseDirectory,i+'_average_expression_matrix.txt')
    average_expression_df.to_csv(avg_file, sep='\t')
```

创建rolypoly前置文件

```R

#install.packages("/share/pub/dengcy/software/rolypoly_0.1.0.tar.gz",repos=NULL,type="source")
library(scPagwas)
library(foreach)
library(rolypoly)
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
Args <- commandArgs(T)
num = print(Args[1])
num<-as.numeric(num)
#for(i in 1:15){
  gwas<-gwass[num]
  cancer<-cancers[num]
  print(gwas)

out_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/1.Gwas_input/{gwas}_gwas_pagwas.RData"))
merge_scexpr<-read.table(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/",cancer,"_average_expression_matrix.txt"))
load(out_file)
ld_path <- "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/rolypoly/LD"
gwas_data<-Pagwas$gwas_data
block_annotation<-block_annotation[!duplicated(block_annotation$label),]
block_annotation$chrom<-sub("chr", "", block_annotation$chrom)
gwas_data$chrom <- sub("chr", "", gwas_data$chrom)
gwas_data$chrom <- as.numeric(gwas_data$chrom)
gwas_data<-gwas_data[,c"chrom","pos","rsid","beta","se","maf")]
block_annotation$chrom <- as.numeric(block_annotation$chrom)
rolypoly_result <- rolypoly_roll(
  gwas_data = gwas_data,
  block_annotation = block_annotation,
  block_data = merge_scexpr,
  ld_folder =ld_path,
  bootstrap_iters = 100
)
result_file= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/rolypoly/{gwas}_rolypoly.RData"))
save(rolypoly_result,file = result_file)


```
服务器计算

```bash
nodes=("in005" "in006" "in007" "in008")
#
#"in005" "in006" "in007" "in008" "in009" "in010"
node_count=${#nodes[@]}
counter=0
# 循环生成.sh文件
for j in {1..15}
do
echo $j
    # 获取当前节点
    node="${nodes[$counter % $node_count]}"
    # 定义.sh文件名
    filename="script_${j}.sh"
    # 创建.sh文件并写入内容
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e roly$j.err
#SBATCH -o roly$j.out
#SBATCH -J roly$j
#SBATCH -w $node
#SBATCH --mem=50000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
Rscript /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/rolypoly/1.r $j
EOF
    # 增加节点计数器
    counter=$((counter + 1))
    sbatch $filename
done
#$node
```

