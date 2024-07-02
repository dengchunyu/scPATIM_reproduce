## magma+scDRS的方法计算的显著性细胞类型

### 输出magma需要的结果
```python
import pandas as pd
import numpy as np
import os
#读取all_cancer_gwas_data.csv文件
gwas_files=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/cancer_gwas_filter.csv")

for i in gwas_files['file_index']:
    print(i)    #输出i
    gwas=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/{i}_isnp.txt".format(i=i),sep=' ')    #读取文件
    #将gwas的pos列转为int类型
    gwas.loc[:,'pos']=gwas.loc[:,'pos'].astype(int)
    magma_Input1=gwas[['rsid','pval']]#提取gwas数据框中的rsid，pval列，赋值给magma_Input1
    #将magma_Input1数据框中加入N列，值为gwas_files中的samplesize列对应的值
    magma_Input1.loc[:,'N']=gwas_files1[gwas_files1['file_index']==i]['samplesize'].values[0]
    magma_Input2=gwas[['rsid','chrom','pos','pos']]#提取gwas数据框中的rsid，chrom，pos，pos列gwas_files中
    magma_Input1.columns=['SNP','P','N']#将magma_Input1的列名改为SNP，P，N
    magma_Input1.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_{i}Input1.txt".format(i=i),sep='\t',index=False)#将magma_Input1的数据框写入文件
    magma_Input2.loc[:,'pos']=magma_Input2.loc[:,'pos'].astype(int)#将magma_Input2的pos列转为int类型
    print(magma_Input2['pos'].dtype)#查看magma_Input2的pos的类型
    magma_Input2.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_{i}Input2.txt".format(i=i),sep='\t',index=False,header=False)#将magma_Input2的数据框写入文件


##输出ldsc需要的结果数据
for i in gwas_files2['gwasdata']:
    gwas=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/{i}_isnp.txt".format(i=i),sep=' ')    #读取文件
    SNP=gwas[['rsid','chrom','pos','REF','ALT','beta','se','p','maf']]
    SNP.columns=['SNP','chrom','pos','A2','A1','beta','se','pval','maf']
    SNP.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/gwas_data/{i}_isnp_ldsc.txt".format(i=i),sep=' ',index=False)

#重新对gwas_files1的filesnames列进行处理，将gwas_files1中的gwasdata列的值赋值给filenames列，对"_isnp_ldsc.txt"文件的A1和A2列进行处理，将A1列的值赋值给A2列，将A2列的值赋值给A1列
for i in gwas_files1['filenames']:
    SNP=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/gwas_data/{i}_isnp_ldsc.txt".format(i=i),sep=' ')
    SNP.columns
    SNP['A2'],SNP['A1']=SNP['A1'],SNP['A2']
    SNP.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/gwas_data/{i}_isnp_ldsc.txt".format(i=i),sep=' ',index=False)

```

所有SNP映射到基因，该窗口从每个基因的上游10kb开始到下游10kb结束**
```shell

magma_filenames=$(cat /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_filenames.txt)
#将magma_filenames.txt文件读取以\t分割，除去第一行，将10列的值赋值给gwass变量
gwass=$(echo "$magma_filenames" | awk -F '\t' 'NR>1{print $10}')
nodes=("in001" "in002" "in003" "in004" "in005" "in006" "in007" "in008")
node_count=${#nodes[@]}
counter=0
#查看文件head
for j in $gwass
do
#查看/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_${j}Input2.txt文件的开头
echo $j
head /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_${j}Input1.txt
done



# 循环生成.sh文件
for j in $gwass
do
    # 获取当前节点
    node="${nodes[$counter % $node_count]}"

    # 定义.sh文件名
    filename="script_${j}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e $j.err
#SBATCH -o $j.out
#SBATCH -J $j
#SBATCH -w $node
#SBATCH --mem=30000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
#对magma_filenames.txt文件中的filenames列进行循环，将每一行的filenames列的值作为文件名进行输入处理
cd /share/pub/dengcy/software/magma/
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_${j}Input2.txt \
--gene-loc /share/pub/dengcy/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/${j}annotated_10kbup_10_down
./magma --bfile /share/pub/dengcy/NCBI/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/magma_${j}Input1.txt ncol=3 \
--gene-annot /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/predata/${j}annotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/result/${j}annotated_10kbup_10down
EOF
    # 提交任务
    sbatch "$filename"

    # 计数器+1
    counter=$((counter + 1))
done
```

```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(org.Hs.eg.db)
library(dplyr)
g2s=toTable(org.Hs.egSYMBOL)
#library(biomaRt)
#读取所有gwas数据的magma结果数据框
gwass<-c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_bladder_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
load("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/g2s.RData")
#df_list<-list()
for(i in gwass){
  df<-read.table(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/result/",i,"annotated_genes.out.csv"),header=T,sep=",")

#将magma_df中的GENE列的名字改为gene_id
colnames(df)[1]<-"gene_id"
magma_df=merge(df,g2s,by="gene_id",all.x=T)
#将magma_df中的gene_id列删除
magma_df<-magma_df[,-1]
colnames(magma_df)[length(colnames(magma_df))]<-"gene_symbol"
magma_df$logp<- -log10(magma_df$P)
write.csv(magma_df,paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/result/",i,"annotated_symbol.out.csv"), quote = F,row.names=F)
}

```

### 计算scDRS细胞类型结果

```bash
source activate mypy

cancers=('HeadandNeck' 'BreastCancer' 'ProstateCancer' 'Melanoma' 'LungCancer' 'BladderCancer' 'KidneyCancer' 'ThyroidCancer' 'OvarianCancer' 'Pancreatic' 'GastricCancer' 'LiverCancer' 'ColorectalCancer' 'EsophagealCancer' 'EndometrialCarcinoma')
gwass=('ieu-b-4912' 'finngen_r7_c3_breast_exallc' 'finngen_r7_c3_prostate_exallc' 'finngen_r7_c3_melanoma' 'finngen_r7_c3_lung_nonsmall_exallc' 'finngen_r7_c3_bladder_exallc' 'finngen_r7_c3_kidney_notrenalpelvis_exallc' 'finngen_r7_c3_thyroid_gland_exallc' 'finngen_r7_c3_ovary_exallc' 'finngen_r7_c3_pancreas_exallc' 'finngen_r7_c3_stomach_exallc' 'bbj-a-158' 'finngen_r7_c3_colorectal' 'bbj-a-117' 'ebi-a-GCST006464')

nodes=("in002" "in003" "in004")
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
#SBATCH -e drs_${cancer}.err
#SBATCH -o drs_${cancer}.out
#SBATCH -J drs_${cancer}
#SBATCH -w ${node}
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
python /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/${cancer}_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime1.0/5.magma/result/${gwas}annotated_symbol.out.csv --top_gene_num 1000 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/scDRS/${gwas}_scDRS_score.csv --group merge_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/scDRS/${gwas}_mergecelltypeP.csv
EOF
  # 提交任务脚本
  sbatch "$filename"
  counter=$((counter + 1))
done

```