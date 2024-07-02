
```python
def un_gz(file_name):
    """ungz zip file"""
    f_name = file_name.replace(".gz", "")
    #获取文件的名称，去掉
    g_file = gzip.GzipFile(file_name)
    #创建gzip对象
    open(f_name, "wb+").write(g_file.read())
    #gzip对象用read()打开后，写入open()建立的文件里。
    g_file.close() #关闭gzip对象


def un_tar(file_name):
    #untar zip file
    tar = tarfile.open(file_name)
    names = tar.getnames()
    if os.path.isdir(file_name + "_files"):
        pass
    else:
        os.mkdir(file_name + "_files")
    #由于解压后是许多文件，预先建立同名文件夹
    for name in names:
        tar.extract(name, file_name + "_files/")
    tar.close()

import scipy
from scipy.sparse import *
import anndata as ad
######################
#colorectalcancer
#######################
###########
#1.E-MTAB-8107
##########
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/E-MTAB-8107/files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/E-MTAB-8107/files')
adatas=[]

for i in files:
    adata = sc.read_csv(i)
    adata = adata.T
    strlist = i.split('.')
    adata.obs['patients']= strlist[0]
    adata.obs['samples']= strlist[0]
    adata.var_names_make_unique()
    adatas.append(adata)

var_names = adatas[0].var_names.union(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.union(adatas[j].var_names)


adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=files)
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/E-MTAB-8107/adata.h5ad"
adata.write(results_file)

###########
#2.GSE164522
##########
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522')

for i in ["GSE164522_CRLM_LN_expression.csv.gz","GSE164522_CRLM_PT_expression.csv.gz","GSE164522_CRLM_metadata.csv.gz"]:
	un_gz(i)
#	
adatas=[]
for i in ["GSE164522_CRLM_LN_expression.csv","GSE164522_CRLM_PT_expression.csv"]:
    adata = sc.read_csv(i)
    adata=adata.T
    strlist = i.split('_')
    adata.obs['patients']= strlist[0]
    adata.obs['samples']= strlist[0] + strlist[2]
    adata.var_names_make_unique()
    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=["GSE164522_CRLM_LN_expression.csv","GSE164522_CRLM_PT_expression.csv"])

adata.obs=adata.obs.drop(labels=['batch'],axis=1)  
del adatas
#############
meta<-read.csv("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/GSE164522_CRLM_metadata.csv")
meta$X<-str_replace_all(meta$X,"-",".")
write.csv(meta,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/GSE164522_CRLM_metadata.csv")
#############
meta=pd.read_csv('GSE164522_CRLM_metadata.csv')
meta.index = meta.iloc[:,0]
mt=[]
for x in meta.index:
    x=x.replace('-', '.')
    mt.append(x)

meta.index=mt
obs = meta.loc[adata.obs_names]
adata.obs['patients']=obs['patient']
adata.obs['celltype_major']=obs['celltype_major']
#3）输出
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/adata.h5ad"
adata.write(results_file)

################
##3.GSE144735
#############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE144735')
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE144735')

for i in files:
	un_gz(i)

adata = sc.read_csv('GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt',delimiter='\t')
meta=pd.read_csv('GSE144735_processed_KUL3_CRC_10X_annotation.txt',sep='\t')
adata = adata.T
# Append metadata
print('Appending metadata')
#2）整合meta数据与矩阵
meta.index=meta.Index
adata.obs = meta

adata.obs['patients']=adata.obs['Patient']
adata.obs['samples']=adata.obs['Sample']
adata.obs=adata.obs.drop(labels=['Index', 'Patient', 'Class', 'Sample', 'Cell_subtype'],axis=1)

#3）输出
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE144735/adata.h5ad"
adata.write(results_file)

################
##4. GSE132257
#############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132257')
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132257')

for i in files:
	un_gz(i)

adata = sc.read_csv('GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt',delimiter='\t')
meta=pd.read_csv('GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt',sep='\t')
adata = adata.T
adata.obs = meta

adata.obs['patients']=adata.obs['Patient']
adata.obs['samples']=adata.obs['Sample']
adata.obs=adata.obs.drop(labels=['Index', 'Patient', 'Class', 'Status', 'Sample'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132257/adata.h5ad"
adata.write(results_file)

############
##5.GSE132465
############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132465')
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132465')

for i in files:
	un_gz(i)

adata = sc.read_csv('GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt',delimiter='\t')
meta=pd.read_csv('GSE132465_GEO_processed_CRC_10X_cell_annotation.txt',sep='\t')
adata = adata.T
adata.obs = meta
adata.obs['patients']=adata.obs['Patient']
adata.obs['samples']=adata.obs['Sample']
adata.obs=adata.obs.drop(labels=['Index', 'Patient', 'Class', 'Sample', 'Cell_subtype'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132465/adata.h5ad"
adata.write(results_file)
##########
#二、整合
##########
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')

adata1 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/E-MTAB-8107/adata.h5ad')
sc.tl.pca(adata1)
sc.external.pp.bbknn(adata1, batch_key='patients')

adata2 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/adata.h5ad')
adata2 = adata2[adata2.obs.celltype_major != 'CD45-']
sc.tl.pca(adata2)
sc.external.pp.bbknn(adata2, batch_key='patients')

adata3 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE144735/adata.h5ad')
sc.tl.pca(adata3)
sc.external.pp.bbknn(adata3, batch_key='patients')

adata4 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132257/adata.h5ad')
sc.tl.pca(adata4)
sc.external.pp.bbknn(adata4, batch_key='patients')

adata5 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132465/adata.h5ad')
sc.tl.pca(adata5)
sc.external.pp.bbknn(adata5, batch_key='patients')

adata_concat = adata1.concatenate(adata2,adata3,adata4,adata5, batch_categories=['GSE164522','E-MTAB-8107', 'GSE144735','GSE132257','GSE132465']) #合并数据集
import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
adata_concat.obs['datasets']=adata_concat.obs.batch

color_map='CMRmap'
cancer='crc'
n_comps=50
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')


sc.pp.pca(adata4)
sc.pp.neighbors(adata4)
sc.tl.umap(adata4)

var_names = adata1.var_names.intersection(adata3.var_names)
var_names = var_names.intersection(adata2.var_names)
var_names = var_names.intersection(adata4.var_names)
var_names = var_names.intersection(adata5.var_names)
adata1=adata1[:,var_names]
adata3=adata3[:,var_names]
adata4=adata4[:,var_names]
adata5=adata5[:,var_names]
sc.tl.ingest(adata1, adata4, obs='Cell_type')
sc.tl.ingest(adata3, adata4, obs='Cell_type')
sc.tl.ingest(adata5, adata4, obs='Cell_type')

adata_concat = adata1.concatenate(adata3,adata4,adata5, batch_categories=['E-MTAB-8107', 'GSE144735','GSE132257','GSE132465']) #合并数据集
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'CRC_integrate_cell.h5ad'

sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch') 
adata_concat.obs['datasets']=adata_concat.obs.batch

adata_concat.write(results_file)

adata_concat = adata_concat[adata_concat.obs.Cell_type != 'Epithelial cells']
adata_concat = adata_concat[adata_concat.obs.Cell_type != 'Stromal cells']


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'CRC_Immune_cell.h5ad'
adata_concat.write(results_file)


sc.tl.pca(adata_concat, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='crc'
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_'+ 'comps_PCA')
#adata = sc.read_h5ad(results_file)
#sc.tl.pca(adata_concat)
#sc.external.pp.bbknn(adata_concat, batch_key='batch')  
#sc.pp.normalize_total(adata_concat, target_sum=1e4) ##标准化
#sc.pp.log1p(adata_concat)
#adata_concat.write(results_file)
```

## 单独处理数据



```python
import scanpy as sc
import pandas as pd
import numpy as np
import os
import gc
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import celltypist
from celltypist import models

#BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/'
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='ColorectalCancer_Immune_cell.h5ad'

adata = sc.read_h5ad(results_file)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig('total_counts_distplot.pdf',format="pdf")
plt.clf()

# sc.pl.violin(adata, 'total_counts')
p2 = sc.pl.violin(adata, "pct_counts_mt")
plt.savefig('pct_counts_mt_violin.pdf',format="pdf")
plt.clf()

p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig('total_counts_genes_scatter.pdf',format="pdf")
plt.clf()

from scipy.stats import median_abs_deviation

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
adata.obs["outlier"] = (is_outlier(adata, "log1p_total_counts", 5) | is_outlier(adata, "log1p_n_genes_by_counts", 5) | is_outlier(adata, "pct_counts_in_top_20_genes", 5))

adata.obs.outlier.value_counts()
#False    1156224
#True       21905
adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)

adata.obs.mt_outlier.value_counts()
print(f"Total number of cells: {adata.n_obs}")
#Total number of cells: 193699
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
#1156224

p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig('total_counts_genes_scatter2.pdf',format="pdf")
plt.clf()

adata.layers["counts"] = adata.X

sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")
adata.write(results_file)
#adata = sc.read_h5ad(results_file)
```

### normalization

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)

sc.pl.umap(adata,color=["total_counts", "pct_counts_mt"],save='_counts')

import seaborn as sns
import matplotlib.pyplot as plt

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata.write(results_file)
```

### annotation

```python
#results_file=os.path.join(BaseDirectory, 'adata_pp_immune.h5ad')
adata = sc.read_h5ad(results_file)
import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")
#model_high = models.Model.load(model="Immune_All_High.pkl")
#predictions_high = celltypist.annotate(adata, model=model_high, majority_voting=True)
#predictions_high_adata = predictions_high.to_adata()
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

predictions_low = celltypist.annotate(adata_celltypist, model=model_low, majority_voting=True)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "majority_voting"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()

sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
# 使用pivot_table计算每个类别中不同细胞类型的计数
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

# 根据每个类别中占比最大的细胞类型来确定细胞名
celltype_names = pivot_table.idxmax(axis=1)

adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)

adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_crc_celltypist_cell_label_fine')
adata = sc.read_h5ad(results_file)
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Endothelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Double-positive thymocytes"]
adata = adata[adata.obs["louvain_1_anno"] != "Fibroblasts"]

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
celltype_names = pivot_table.idxmax(axis=1)
adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_crc_celltypist_cell_label_fine')

adata.write(results_file)

cancer='ColorectalCancer'
adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/"+ cancer + "_Immune_cell.h5ad")
adata=adata[adata.obs.first_celltype_annotation != 'NK']
sc.pl.umap(adata,color=["first_celltype_annotation"],save='_ColorectalCancer_first_celltype_annotation')
adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].cat.remove_unused_categories()

#marker_dotplot(cancer="ColorectalCancer")
results_file='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/ColorectalCancer_Immune_cell.h5ad'

adata.write(results_file)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
adata = sc.read_h5ad(results_file)

```