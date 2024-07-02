# scPagwas原始结果的细胞类型比较

```R
library(Seurat)
library(Matrix)
library(scPagwas)
Args <- commandArgs(T)
num = print(Args[1])
#'HeadandNeck',
cancers=c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass=c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.2.scPagwas_orignal")
Args <- commandArgs(T)
num = print(Args[1])
num<-as.numeric(num)
gwas<-gwass[num]
cancer<-cancers[num]
print(gwas)
setwd('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.2.scPagwas_orignal')
scdata_file=paste0('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/',i,'_scdata.rds')
scdata<-readRDS(scdata_file)
gwasfile= file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/{gwas}_isnp.txt"))
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =gwasfile,
                     Single_data =scdata_file,
                     output.prefix="",
                     output.dirs=cancer,
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     block_annotation = block_annotation,
                     iters_singlecell = 0,
                     chrom_ld = chrom_ld,
                     seurat_return=F,
                     singlecell=F, # Whether to run the singlecell process.
                     celltype=T)
```

```bash
nodes=("in006" "in007" "in008")
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
#SBATCH -e scPagwas_$j.err
#SBATCH -o scPagwas_$j.out
#SBATCH -J scPagwas_$j
#SBATCH -w $node
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
Rscript /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.2.scPagwas_orignal/1.r $j
EOF
    # 增加节点计数器
    counter=$((counter + 1))
    sbatch $filename
done

```
