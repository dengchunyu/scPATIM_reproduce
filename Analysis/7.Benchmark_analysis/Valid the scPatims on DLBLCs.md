# 基于DLBLCs数据进行scPagwas-eQTLs框架的验证

## DLBLC 全部gwas数据处理

```r
#finngen_R9_C3_DLBCL_EXALLC
GWAS_raw <-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/rawdata/finngen_R9_C3_DLBCL_EXALLC.gz")
colnames(GWAS_raw )<-c("chrom","pos","REF","ALT","rsid","nearest_genes","p","lp","beta","se","maf","maf_case","maf_control")
GWAS_raw$N<-287137
write.table(GWAS_raw,file= "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_gwas_data.txt",row.names=F,quote=F,sep="\t")
gwas_data<-GWAS_raw[!GWAS_raw$rsid %in% GWAS_raw$rsid[duplicated(GWAS_raw$rsid)],]
write.table(gwas_data,file= "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_gwas_data.txt",row.names=F,quote=F,sep="\t")


magma1 <- gwas_data[,c('rsid','chrom','pos','pos')]
magma2 <- gwas_data[,c('rsid','p','N')]
colnames(magma2)<-c('SNP','P','N')
write.table(magma1,"magma1.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(magma2,"magma2.txt",sep="\t",quote=F,row.names=F,col.names=F)
```

开始计算LD

```bash
#提取snp
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/
awk  '{print $5 }' DLBCL_gwas_data.txt  > SNP_list.txt
mkdir tempfile

for i in $(seq 1 22)  
do
/share/pub/dengcy/software/plink-1.07-x86_64/plink --bfile /share/pub/dengcy/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i --extract SNP_list.txt --noweb --make-bed --out ./tempfile/1000G.EUR.QC.${i}_filtered
done

#计算LD信息
for i in $(seq 1 22)
do 
echo $i
/share/pub/dengcy/software/plink-1.07-x86_64/plink --bfile ./tempfile/1000G.EUR.QC.${i}_filtered --indep-pairwise 50 5 0.8 --out  ./tempfile/${i}_plink_LD0.8
done

cat ./tempfile/*.prune.in > ./tempfile/SNP_LD0.8.prune.txt
```

### 

```r
library(bigreadr)
gwas_data <- bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result//DLBCL_gwas_data.txt")
snps<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/tempfile/SNP_LD0.8.prune.txt")
snps<-unique(snps[,1])
#判断gwas_data$rsid是否有重复，如果有重复，就删除重复的行
gwas_data2<-gwas_data[gwas_data$rsid %in% snps,]
write.table(gwas_data2,file= "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_prune_gwas_data.txt",row.names=F,quote=F,sep="\t")
#eqtls filter
library(readr)
library(dplyr)
immsnp<-read_table('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_filter5e2.tsv',col_names = F)

snps<-intersect(gwas_data$rsid,immsnp$X2)
print(length(snps))
#[1] 2440148
dim(gwas_data)
#[1] 2440148      14
gwas_data3<-gwas_data[gwas_data$rsid %in% snps,]
write.table(gwas_data3,file= "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt",row.names=F,quote=F,sep="\t")

magma1 <- gwas_data[,c('rsid','chrom','pos','pos')]
magma2 <- gwas_data[,c('rsid','p','N')]
colnames(magma2)<-c('SNP','P','N')
#将magma1保存到/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma1/中
write.table(magma1,paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/magma1/DLBCL_allgwas.txt"),sep="\t",quote=F,row.names=F,col.names=F)
#将magma2保存到/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma2/中
write.table(magma2,paste0("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/magma2/DLBCL_allgwas.txt"),sep="\t",quote=F,row.names=F,col.names=F)
```

```bash
cd /share/pub/dengcy/software/magma/
./magma --annotate window=10,10 --snp-loc /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/magma1/DLBCL_allgwas.txt \
--gene-loc /share/pub/dengcy/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded \
--out /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/DLBCL_allgwasannotated_10kbup_10_down
./magma --bfile /share/pub/dengcy/NCBI/g1000_eur/g1000_eur \
--pval /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/magma2/DLBCL_allgwas.txt ncol=3 \
--gene-annot /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/DLBCL_allgwasannotated_10kbup_10_down.genes.annot \
--out /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/magma_result/DLBCL_allgwasannotated_10kbup_10down
```

## PBMC

### 整合pbmc的h5ad文件

```python
import scanpy as sc
import pandas as pd
import numpy as np
import os
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/")
#scanpy读取NM_Healthy_pbmc_count.csv
adata = sc.read_csv("/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc_count.csv")
#转置
adata = adata.transpose()
#读取NM_Healthy_pbmc_meta.csv
meta = pd.read_csv("/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc_meta.csv")
adata.obs = meta
#输出adata到h5ad文件中
###处理
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata)
#注释聚类结果
adata.obs['annotation'] = adata.obs['leiden'].apply(lambda x: adata.obs[adata.obs['leiden']==x]['initial_clustering'].value_counts().index[0])

sc.pl.umap(adata, color=['annotation'],save="DLBCL_annotation",wspace=0.5)
#将annotation==B cells，设置为阳性细胞pos，其余设置为阴性细胞neg
adata.obs['label'] = adata.obs['annotation'].apply(lambda x: 1 if x=="B_cell" else 0)
adata.obs['annotation'].value_counts()
adata.obs['label'].value_counts()

adata.write("NM_Healthy_pbmc.h5ad")

```

### 计算scPagwas
这里为了和pipeline进行比较，将使用kegg进行scPagwas的计算
这里将gwas数据进行eqtls的处理其他的不变。
```R
library(scPagwas)
library(Seurat)
library(dplyr)

setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result")
scdata <- readRDS("/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc.rds")

load(scmat_path)
output.dirs="DLBCL_scPagwas"
gwas_data <- "DLBCL_eqtls_gwas_data.txt"
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =gwas_data,
                     Single_data ="/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc.rds",
                     output.prefix="",
                     output.dirs=output.dirs,
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     block_annotation = block_annotation,
                     iters_singlecell = 0,
                     chrom_ld = chrom_ld,# The LD data is provided by package.
                     singlecell=T, # Whether to run the singlecell process.
                     celltype=T)
```
sbatch scPagwas.sh


### scPatim的计算
```bash
#!/bin/bash
#SBATCH -e DLBCL.err
#SBATCH -o DLBCL.out
#SBATCH -J DLBCL
#SBATCH -w in009
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
cd /share/pub/dengcy/Cancer_Gwas/src3.0/scPagwas_eqtls/
scdata_ad='/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc.rds'
gwas_ad='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt'
setwd='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/'
iters = 0
out='DLBCL_eqtls'
Rscript 1.scPagwas_eqtls_pipeline.r $scdata_ad $gwas_ad $setwd $out $iters
```

下面是1.scPagwas_eqtls_pipeline.r中对于基因的计算：
```R
Random_PCC<- function(gPas,datamat,seed=1234,random_times=100,select_num=10000){
    message("Start to run the corSparse function in random!")
    set.seed(seed)
    sparse_cor_list<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(datamat),select_num)
    gPas_select <- data.matrix(gPas[index])
    sparse_cor <- scPagwas::corSparse(
      X = t(datamat[,index]),
      Y = gPas_select
    )
    colnames(sparse_cor) <- "PCC"
    sparse_cor[is.nan(sparse_cor)] <- 0
    sparse_cor[is.na(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
    }
    sparse_cor_list<-as.data.frame(sparse_cor_list)
    sparse_cor<- apply(sparse_cor_list,1,function(x) mean(x, na.rm = TRUE))
    return(data.matrix(sparse_cor))
}
library(scPagwas)
library(Seurat)
library(dplyr)
setwd('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result')
out='DLBCL_eqtls'
output.dirs<-out
output.prefix<-""
scdata_ad='/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc.rds'
Single_data<-readRDS(scdata_ad)

load(paste0("./",out,"/Pagwas.RData"))
scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
Pagwas$data_mat <- GetAssayData(Single_data, slot = "data", assay = "RNA")
Pagwas$data_mat <- as_matrix(Pagwas$data_mat)


pcc <- Random_PCC(gPas=scPagwas.gPAS.score,datamat=Pagwas$data_mat,seed=1234,random_times=100,select_num=floor(ncol(Pagwas$Pathway_single_results)/2))
colnames(pcc)<-'PCC'

top5_index<-order(scPagwas.gPAS.score,decreasing=T)[1:(0.05*(length(scPagwas.gPAS.score)))]
    top5_cellnames<-names(scPagwas.gPAS.score)[top5_index]
    top5_index2<-names(scPagwas.gPAS.score) %in% top5_cellnames
    p_list<-apply(Pagwas$data_mat,1,function(x){
        result <- wilcox.test( x[top5_index2], x[!top5_index2], alternative = "greater")
        adjusted_p_value <- p.adjust(result$p.value, method = "bonferroni")
        return(adjusted_p_value)
        })
pcc<-as.data.frame(pcc)
pcc$pvalue <- p_list
#矫正pvalue
pcc$adj_pvalue <- p.adjust(pcc$pvalue, method = "bonferroni")
pcc$adj_logp<- -log10(pcc$adj_pvalue)
pcc$adj_logp[!is.finite(pcc$adj_logp)]<- max(pcc$adj_logp[is.finite(pcc$adj_logp)])+1

pcc$weight_pcc<- pcc$adj_logp * pcc$PCC

utils::write.csv(pcc,file = paste0("./", output.dirs, "/",output.prefix,"_gene_PCC.csv"), quote = F)

n_topgenes=500
assay="RNA"
scPagwas_topgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = T)[1:n_topgenes]]
scPagwas_downgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = F)[1:n_topgenes]]

Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))
Single_data<-Single_data[,names(scPagwas.gPAS.score)]

a <- data.frame(
#scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
#scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
scPagwas.gPAS.score = scPagwas.gPAS.score
)
utils::write.csv(a,file = paste0("./", output.dirs, "/",output.prefix,"_singlecell_scPagwas_score_pvalue.Result.csv"),quote = F)
```

### 计算pbmc 不同方法的score

```bash
source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
scRNA_h5ad_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad
magma_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_gene_pvalue.csv
scpagwas_file=/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_scPagwas/_gene_PCC.csv
scpatim_file=/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls/_gene_PCC.csv

python scDRS_pipeline_nocelltype.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad --gene_file $scpagwas_file --top_gene_num 500 --n_ctrl 500  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_pbmc_orignal_scPagwas.txt --weight_pcc PCC

python scDRS_pipeline_nocelltype.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad --gene_file $scpatim_file --top_gene_num 500 --n_ctrl 500  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_pbmc_scPatim.csv --weight_pcc weight_pcc

python scHII_score.py --scRNA_h5ad_file $scRNA_h5ad_file --scHII_gene_file $magma_file --top_gene_num 500 --n_ctrl 500 --scHII_score_file /share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_pbmc/magma_score.txt
```

### pbmc整合结果并且画图

```python

import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result")
def readdata():
  adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")
  adata.obs=adata.obs.drop(labels=['scPagwas_eqtls_norm_score', 'scPagwas_eqtls_mc_pval', 'scPagwas_norm_score', 'scPagwas_mc_pval', 'scDRS_norm_score', 'scDRS_mc_pval',],axis=1)
  scPagwas_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_pbmc_orignal_scPagwas.txt"
  scpatim_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_pbmc_scPatim.csv"
  magma_file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_pbmc/magma_score.txt"
  scPagwas_score = pd.read_csv(scPagwas_file)
  scPagwas_score.index = scPagwas_score.iloc[:, 0]
  scPagwas_score = scPagwas_score.iloc[:, 1:3]
  scPagwas_score.columns = ['scPagwas_' + i for i in scPagwas_score.columns]
  scpatim_score = pd.read_csv(scpatim_file)
  scpatim_score.index = scpatim_score.iloc[:, 0]
  scpatim_score = scpatim_score.iloc[:, 1:3]
  scpatim_score.columns = ['scPatim_' + i for i in scpatim_score.columns]
  scpatim_score.index=adata.obs.index
  scPagwas_score.index=adata.obs.index
  magma_score = pd.read_csv(magma_file)
  magma_score.index = magma_score.iloc[:, 0]
  magma_score = magma_score.iloc[:, 1:3]
  magma_score.columns = ['scDRS_' + i for i in magma_score.columns]
  magma_score.index=adata.obs.index
  adata.obs = pd.concat([adata.obs,scpatim_score, scPagwas_score,magma_score], axis=1)
  return adata

adata=readdata()
roc_auc_score(adata.obs['label'],adata.obs['scPatim_norm_score'])
roc_auc_score(adata.obs['label'],adata.obs['scPagwas_norm_score'])
roc_auc_score(adata.obs['label'],adata.obs['scDRS_norm_score'])

sc.pl.umap(adata, color=['scPagwas_norm_score','scPatim_norm_score','scDRS_norm_score'],save="DLBCL_scPagwas_scPatim_scDRS_score",wspace=0.5)


plt.clf()
def plot_roc(adata):
    roc_auc_score(adata.obs['label'],adata.obs['scPatim_norm_score'])
    roc_auc_score(adata.obs['label'],adata.obs['scPagwas_norm_score'])
    roc_auc_score(adata.obs['label'],adata.obs['scDRS_norm_score'])
    plt.figure(figsize=(6, 6))
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scPatim_norm_score'])
    plt.plot(fpr,tpr,label="scPATIM.TRS.Score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scPatim_norm_score']))[0:5])
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scPagwas_norm_score'])
    plt.plot(fpr,tpr,label="scPagwas.TRS.Score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scPagwas_norm_score']))[0:5])
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scDRS_norm_score'])
    plt.plot(fpr,tpr,label="scDRS_score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scDRS_norm_score']))[0:5])
    plt.legend()
    plt.savefig("pbmc_scPagwas_scPatim_roc.pdf")
    plt.close()

plot_roc(adata)

#输出adata
adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")
#输出obs
adata.obs.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_pbmc/DLBCL_scPtim_scPagwas_magma.csv")
```


## DLBLC疾病单细胞数据计算所有的免疫细胞的得分

### 读取DLBCL单细胞数据

GSE182434
```shell
#解压/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_raw_count_matrix.txt.gz
gunzip /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_raw_count_matrix.txt.gz
#解压/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_cell_annotation.txt.gz
gunzip /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_cell_annotation.txt.gz
```

输出rds文件
```r
#读取GSE182434_raw_count_matrix.txt.gz
library(Seurat)
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/")
#读取/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_raw_count_matrix.txt
data<-read.table("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_raw_count_matrix.txt",sep="\t",header = T,row.names = 1)
metadata<-read.table("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_cell_annotation.txt",sep="\t",header = T)
id<-intersect(colnames(data),metadata$ID)
data<-data[,id]
metadata<-metadata[metadata$ID %in% id,]
rownames(metadata)<-metadata$ID
#构造Seurat对象
DLBCL<-CreateSeuratObject(counts= data,project = "DLBCL",min.cells = 3,min.features = 200,meta.data = metadata)
Idents(DLBCL)<-DLBCL$CellType
#预处理
DLBCL <- NormalizeData(DLBCL, normalization.method = "LogNormalize", scale.factor = 10000)
DLBCL<-FindVariableFeatures(DLBCL,selection.method = "vst",nfeatures = 2000)
DLBCL<-ScaleData(DLBCL)
DLBCL<-RunPCA(DLBCL,features = VariableFeatures(object = DLBCL))
#保存DLBCL到rds文件中
saveRDS(DLBCL,"DLBCL.rds")

DLBCL<-readRDS("DLBCL.rds")
#提取B cells
DLBCL_B_cells<-DLBCL[,DLBCL$CellType=="B cells"]
#聚类处理
#pca
DLBCL_B_cells<-RunPCA(DLBCL_B_cells,features = VariableFeatures(object = DLBCL_B_cells))
DLBCL_B_cells<-FindNeighbors(DLBCL_B_cells, dims = 1:30)
DLBCL_B_cells<-FindClusters(DLBCL_B_cells, resolution = 0.2)
Idents(DLBCL_B_cells)<-DLBCL_B_cells$seurat_clusters
#保存DLBCL_B_cells到rds文件中
saveRDS(DLBCL_B_cells,"DLBCL_B_cells.rds")
DLBCL_B_cells<-readRDS("DLBCL_B_cells.rds")
#去批次
library(harmony)
DLBCL_B_cells<-RunHarmony(object = DLBCL_B_cells,group.by.vars = "Patient")
#重新聚类
DLBCL_B_cells<-FindNeighbors(DLBCL_B_cells, dims = 1:30)
DLBCL_B_cells<-FindClusters(DLBCL_B_cells, resolution = 0.2)
#umap
DLBCL_B_cells<-RunUMAP(DLBCL_B_cells, dims = 1:30)
Idents(DLBCL_B_cells)<-DLBCL_B_cells$seurat_clusters
#可视化
pdf("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B/DLBCL_B_cells_seurat_clusters.pdf")
DimPlot(DLBCL_B_cells, reduction = "umap",label=T,group.by = "Patient")
dev.off()
```


输出h5ad文件
```python
import scanpy as sc
import pandas as pd
import numpy as np
import os
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/")
#读取GSE182434_raw_count_matrix.txt.gz
adata = sc.read_text("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_raw_count_matrix.txt.gz")
#转置
adata = adata.transpose()
#读取注释文件GSE182434_cell_annotation.txt.gz
meta = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/DLBLC/GSE182434_cell_annotation.txt.gz",sep="\t")
adata.obs = meta
#选择Tissue=="DLBCL"的细胞
#adata = adata[adata.obs['Tissue']=="DLBCL",:]
#查看基因
adata.var
#输出adata到h5ad文件中
adata.write("DLBCL.h5ad")
adata= sc.read_h5ad("DLBCL.h5ad")
#提取B cells
adata.obs['CellType'].value_counts()
adata = adata[adata.obs['CellType']=="B cells",:]
#保存adata到h5ad文件中
adata.write("DLBCL_B_cells.h5ad")
#标准化
sc.pp.normalize_total(adata, target_sum=1e4)
#对adata进行log1p
sc.pp.log1p(adata)
#根据Patient进行去批次处理
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata)
#去批次
sc.pp.combat(adata, key='Patient')
#重新计算umap
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata)
#保存adata到h5ad文件中
adata.write("DLBCL_B_cells.h5ad")

adata.obs.Tissue.value_counts()
#FL        4816
#DLBCL     3436
#Tonsil    1663
#新建一列cell_type,paste Tissue 和TumorNormal
adata.obs['Tissue'] = adata.obs['Tissue'].astype(str)
adata.obs['TumorNormal'] = adata.obs['TumorNormal'].astype(str)
adata.obs['cell_type'] = adata.obs['Tissue']+"_"+adata.obs['TumorNormal']
adata.obs['cell_type'].value_counts()
#可视化，TumorNormal
adata.obs.Tissue = adata.obs.Tissue.astype(str)
adata.obs.TumorNormal = adata.obs.TumorNormal.astype(str)
adata.obs['cell_type'] = adata.obs.Tissue + "_" + adata.obs.TumorNormal

sc.pl.umap(adata, color=['Patient','cell_type'],save="DLBCL_B_cell_type",wspace=0.5)
adata.obs.cell_type.value_counts()
adata.obs.cell_type = adata.obs.cell_type.astype(str)
#将cell_type==DLBCL_Tumor，设置为阳性细胞pos，其余设置为阴性细胞neg
adata.obs['label'] = adata.obs.cell_type.apply(lambda x: 1 if x=="DLBCL_Tumor" else 0)
adata.write_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad")

```

### scPagwas计算

```R
library(scPagwas)
library(Seurat)
library(dplyr)
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result")

output.dirs="DLBCL_B_original_scPagwas"
gwas_data <- "DLBCL_eqtls_gwas_data.txt"
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =gwas_data,
                     Single_data ="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds",
                     output.prefix="",
                     output.dirs=output.dirs,
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     block_annotation = block_annotation,
                     iters_singlecell = 0,
                     chrom_ld = chrom_ld,
                     seurat_return=F,
                     singlecell=T, # Whether to run the singlecell process.
                     celltype=T)
save(Pagwas,file=paste0("./",output.dirs,"/Pagwas.RData"))

load(paste0("./",output.dirs,"/Pagwas.RData"))
output.prefix<-""
scdata_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
Single_data<-readRDS(scdata_ad)

#load(paste0("./",out,"/Pagwas.RData"))
scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
Pagwas$data_mat <- GetAssayData(Single_data, slot = "data", assay = "RNA")
Pagwas$data_mat <- as_matrix(Pagwas$data_mat)
pcc <- Random_PCC(gPas=scPagwas.gPAS.score,datamat=Pagwas$data_mat,seed=1234,random_times=100,select_num=floor(dim(Pagwas$data_mat)[2]/2))
colnames(pcc)<-'PCC'

top5_index<-order(scPagwas.gPAS.score,decreasing=T)[1:(0.05*(length(scPagwas.gPAS.score)))]
top5_cellnames<-names(scPagwas.gPAS.score)[top5_index]
top5_index2<-names(scPagwas.gPAS.score) %in% top5_cellnames
p_list<-apply(Pagwas$data_mat,1,function(x){
    a1<- na.omit(x[top5_index2])
    a2<- na.omit(x[!top5_index2])
    result <- wilcox.test(a1 ,a2 , alternative = "greater")
    adjusted_p_value <- p.adjust(result$p.value, method = "bonferroni")
    return(adjusted_p_value)
    })
pcc<-as.data.frame(pcc)
pcc$PCC[is.nan(pcc$PCC)]<-0
pcc$pvalue <- p_list
#矫正pvalue
pcc$adj_pvalue <- p.adjust(pcc$pvalue, method = "bonferroni")
pcc$adj_logp<- -log10(pcc$adj_pvalue)
pcc$adj_logp[!is.finite(pcc$adj_logp)]<- max(pcc$adj_logp[is.finite(pcc$adj_logp)])+1
#pcc$logp<- -log10(pcc$pvalue)
pcc$weight_pcc<- pcc$adj_logp * pcc$PCC
utils::write.csv(pcc,file = paste0("./", output.dirs, "/",output.prefix,"_gene_PCC.csv"), quote = F)
```

sbatch scPagwas_B.sh

### scpatim

```bash
####
#!/bin/bash
#SBATCH -e DLBCL3.err
#SBATCH -o DLBCL3.out
#SBATCH -J DLBCL3
#SBATCH -w in008
#SBATCH --mem=150000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
cd /share/pub/dengcy/Cancer_Gwas/src3.0/scPagwas_eqtls/
scdata_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
gwas_ad='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt'
setwd='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/'
iters = 0
out='DLBCL_B_scPATIM'
Rscript 1.scPagwas_eqtls_pipeline.r $scdata_ad $gwas_ad $setwd $out $iters
```

```R
Random_PCC<- function(gPas,datamat,seed=1234,random_times=100,select_num=10000){
    message("Start to run the corSparse function in random!")
    set.seed(seed)
    sparse_cor_list<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(datamat),select_num)
    gPas_select <- data.matrix(gPas[index])
    sparse_cor <- scPagwas::corSparse(
      X = t(datamat[,index]),
      Y = gPas_select
    )
    colnames(sparse_cor) <- "PCC"
    #sparse_cor[is.nan(sparse_cor)] <- 0
    #sparse_cor[is.na(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
    }
    sparse_cor_list<-as.data.frame(sparse_cor_list)
    sparse_cor<- apply(sparse_cor_list,1,function(x) mean(x, na.rm = TRUE))
    return(data.matrix(sparse_cor))
}
library(scPagwas)
library(Seurat)
library(dplyr)
setwd('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/')
out='DLBCL_B_scPATIM'
output.dirs<-out
output.prefix<-""
scdata_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
Single_data<-readRDS(scdata_ad)

load(paste0("./",out,"/Pagwas.RData"))
scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
names(scPagwas.gPAS.score)<-colnames(Single_data)
Pagwas$data_mat <- GetAssayData(Single_data, slot = "data", assay = "RNA")
Pagwas$data_mat <- as_matrix(Pagwas$data_mat)
pcc <- Random_PCC(gPas=scPagwas.gPAS.score,datamat=Pagwas$data_mat,seed=1234,random_times=100,select_num=floor(ncol(Pagwas$Pathway_single_results)/2))
colnames(pcc)<-'PCC'

top5_index<-order(scPagwas.gPAS.score,decreasing=T)[1:(0.05*(length(scPagwas.gPAS.score)))]
top5_cellnames<-colnames(Single_data)[top5_index]
top5_index2<-names(scPagwas.gPAS.score) %in% top5_cellnames
p_list<-apply(Pagwas$data_mat,1,function(x){
    a1<- na.omit(x[top5_index2])
    a2<- na.omit(x[!top5_index2])
    result <- wilcox.test(a1 ,a2 , alternative = "greater")
    adjusted_p_value <- p.adjust(result$p.value, method = "bonferroni")
    return(adjusted_p_value)
    })

pcc<-as.data.frame(pcc)
pcc$PCC[is.nan(pcc$PCC)]<-0
pcc$pvalue <- p_list
#矫正pvalue
pcc$adj_pvalue <- p.adjust(pcc$pvalue, method = "bonferroni")
pcc$adj_logp<- -log10(pcc$adj_pvalue)
pcc$adj_logp[!is.finite(pcc$adj_logp)]<- max(pcc$adj_logp[is.finite(pcc$adj_logp)])+1
pcc$weight_pcc<- pcc$adj_logp * pcc$PCC

utils::write.csv(pcc,file = paste0("./", output.dirs, "/",output.prefix,"_gene_PCC.csv"), quote = F)

n_topgenes=500
assay="RNA"
scPagwas_topgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = T)[1:n_topgenes]]
scPagwas_downgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = F)[1:n_topgenes]]

Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))
Single_data<-Single_data[,names(scPagwas.gPAS.score)]

a <- data.frame(
scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
scPagwas.gPAS.score = scPagwas.gPAS.score
)
utils::write.csv(a,file = paste0("./", output.dirs, "/",output.prefix,"_singlecell_scPagwas_score_pvalue.Result.csv"),quote = F)
```

### scPATIM

```bash
#!/bin/bash
#SBATCH -e DLBCL.err
#SBATCH -o DLBCL.out
#SBATCH -J DLBCL
#SBATCH -w in009
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate R4.2
cd /share/pub/dengcy/Cancer_Gwas/src/scHII_pipeline
scdata_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.rds'
gwas_ad='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/gwasdata/DLBCL_gwas_data.txt'
setwd='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/'
iters = 100
out='DLBCL_B_scPagwas'
Rscript 1.scPagwas_scRNAseq.r $scdata_ad $gwas_ad $setwd $out $iters
```

## 基于结果基因计算得分

```bash
conda deactivate
source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
scRNA_h5ad_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad

### scpagwas
scPagwas_file=/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_original_scPagwas/_gene_PCC.csv
scpagwas_score_file=/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_B_pbmc_original_scPagwas.txt

python scDRS_pipeline_nocelltype.py  --scRNA_h5ad_file $scRNA_h5ad_file --gene_file $scPagwas_file --top_gene_num 200 --n_ctrl 200 --score_file $scpagwas_score_file --weight_pcc weight_pcc

###scpatim
cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
scPagwas_eqtls_file=/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/_gene_PCC.csv

python scDRS_pipeline_nocelltype.py --scRNA_h5ad_file $scRNA_h5ad_file --gene_file $scPagwas_eqtls_file --top_gene_num 200 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/DBCLC_B_scpatim.csv --weight_pcc weight_pcc

#3.magma_score
cd /share/pub/dengcy/Cancer_Gwas/src1.0/scHII_pipeline
out='DLBCL_B'
setwd='/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/'
#第一行加gene
magma_file="${setwd}magma_result/DLBCL_allgwasannotated_10kbup_10down.genes.out"
#只提取magma_file的第1,9列并输出到DLBCL_allgwasannotated_10kbup_10down.genes.out.csv
awk -F" " '{print $1","$9}' ${magma_file} > ${magma_file}.csv
magma_file2="${setwd}/magma_result/DLBCL_allgwasannotated_genesymbol.out.csv"
source activate R4.2
Rscript magma_idconvert.r $magma_file.csv $magma_file2
magma_file="${setwd}/magma_result/DLBCL_allgwasannotated_genesymbol.out.csv"
scHII_gene_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_gene_pvalue.csv
scHII_score_file=/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B/magma_score.txt
python scHII_score.py --scRNA_h5ad_file $scRNA_h5ad_file --scHII_gene_file $scHII_gene_file --top_gene_num 500 --n_ctrl 500 --scHII_score_file $scHII_score_file
```

### 可视化结果

```python
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result")

#########查看gpas得分密度
gPAS_score = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_scPagwas/_singlecell_scPagwas_score_pvalue.Result.csv",sep=",")
gPAS_score.index = gPAS_score['Unnamed: 0']
#删除第一列
gPAS_score = gPAS_score.drop(['Unnamed: 0'],axis=1)
gPAS_score = gPAS_score['scPagwas.gPAS.score']
gPas_density(gPas_score=gPAS_score,filen='gPAS_score_density.pdf')


adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad")
adata.obs.index = adata.obs['ID']
adata.obs = adata.obs.drop(columns=['scHII_mc_pval', 'scHII_norm_score', 'scHII_down_mc_pval', 'scHII_down_norm_score', 'scPagwas_down_mc_pval', 'scPagwas_down_norm_score', 'scPagwas_mc_pval', 'scPagwas_norm_score', 'magma_mc_pval', 'magma_norm_score', 'scPatim_norm_score', 'scPatim_mc_pval', 'scPagwas_norm_score', 'scPagwas_mc_pval', 'scDRS_norm_score', 'scDRS_mc_pval'])
#scPagwas_score

def readdata():
  scPagwas_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DBCLC_B_pbmc_original_scPagwas.txt"
  scpatim_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/scDRS_score_top_500_genes.csv"
  magma_file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B/magma_score.txt"
  scPagwas_score = pd.read_csv(scPagwas_file)
  scPagwas_score.index = scPagwas_score.iloc[:, 0]
  scPagwas_score = scPagwas_score.iloc[:, 1:3]
  scPagwas_score.columns = ['scPagwas_' + i for i in scPagwas_score.columns]
  scpatim_score = pd.read_csv(scpatim_file)
  scpatim_score.index = scpatim_score.iloc[:, 0]
  scpatim_score = scpatim_score.iloc[:, 1:3]
  scpatim_score.columns = ['scPatim_' + i for i in scpatim_score.columns]
  scpatim_score.index=adata.obs.index
  scPagwas_score.index=adata.obs.index
  magma_score = pd.read_csv(magma_file)
  magma_score.index = magma_score.iloc[:, 0]
  magma_score = magma_score.iloc[:, 1:3]
  magma_score.columns = ['scDRS_' + i for i in magma_score.columns]
  magma_score.index=adata.obs.index
  adata.obs = pd.concat([adata.obs,scpatim_score, scPagwas_score,magma_score], axis=1)
  return adata

adata = readdata()

roc_auc_score(adata.obs['label'],adata.obs['scPatim_norm_score'])
roc_auc_score(adata.obs['label'],adata.obs['scPagwas_norm_score'])
roc_auc_score(adata.obs['label'],adata.obs['scDRS_norm_score'])

plt.clf()
def plot_roc(adata):
    plt.figure(figsize=(6, 6))
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scPatim_norm_score'])
    plt.plot(fpr,tpr,label="scpatim.TRS.Score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scPatim_norm_score']))[0:5])
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scPagwas_norm_score'])
    plt.plot(fpr,tpr,label="scPagwas.TRS.Score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scPagwas_norm_score']))[0:5])
    fpr,tpr,thresholds = roc_curve(adata.obs['label'],adata.obs['scDRS_norm_score'])
    plt.plot(fpr,tpr,label="scDRS_score,auc="+str(roc_auc_score(adata.obs['label'],adata.obs['scDRS_norm_score']))[0:5])
    plt.legend()
    plt.savefig("DLBLC_B_scPagwas_scPatim_roc.pdf")
    plt.close()
plot_roc(adata)

sc.pl.umap(adata, color=['scPagwas_norm_score','scPatim_norm_score','scDRS_norm_score'],save="DLBCL_B_scPagwas_scPatim_scDRS_score",wspace=0.5)

#输出adata
adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")
#输出obs
adata.obs.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_pbmc/DLBCL_scpatim_scPagwas_magma.csv")

def gPas_density(gPas_score,filen):
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.clf()
    plt.figure(figsize=(8, 8))
    sns.kdeplot(gPas_score, fill=True)
    mean = np.mean(gPas_score)
    plt.axvline(mean, linestyle='--', color='red')
    plt.xlabel('gPAS score')
    plt.ylabel('Density')
    plt.title('Density plot of gPAS score')
    plt.savefig(filen)
    plt.close()
gPAS_score = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_prune_gwas/_singlecell_scPagwas_score_pvalue.Result.csv",sep=",")
gPAS_score.index = gPAS_score['Unnamed: 0']
#删除第一列
gPAS_score = gPAS_score.drop(['Unnamed: 0'],axis=1)
gPAS_score = gPAS_score['scPagwas.gPAS.score']
gPas_density(gPas_score=gPAS_score,filen='DLBLC_B_gPAS_score_density.pdf')
```
