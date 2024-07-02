```python
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer
def re_size_adata(adata,var_names):
    var_name=var_names.difference(adata.var_names)
    len(var_name)
    n1=len(var_name)
    n2=len(adata.obs_names)
    z=np.zeros((n2,n1))
    z=csc_matrix(z,dtype='float64')
    if isspmatrix_csc(adata.X):
        X = scipy.sparse.hstack((adata.X,z))
    else:
        X = csc_matrix(adata.X,dtype='float64')
        X = hstack((X,z))
    var_names= adata.var_names.union(var_name)
    var = pd.DataFrame(index=var_names)
    adata=ad.AnnData(X,obs=adata.obs,var=var,dtype='float64')
    return adata

import anndata as ad
######################
#GSE146026 只有上皮细胞
############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/GSE146026')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/GSE146026')
for x in files:
	un_gz(x)

#adata=pd.read_csv('GSE146026_Izar_HGSOC_ascites_10x_log.tsv',sep='\t')

adatas = []
'GSE146026_Izar_HGSOC_ascites_10x_log.tsv'

X = pd.read_csv('GSE146026_Izar_HGSOC_ascites_10x_log.tsv',sep='\t',index_col=0,header=1)
n_obs=len(X.columns)
n_vars = len(X.index)
obs = pd.DataFrame()
obs['sample'] = np.tile(i,n_obs)
obs.index=X.columns
var_names=X.index
var = pd.DataFrame(index=var_names)
X=X.T
adata=ad.AnnData(X,obs=obs,var=var,dtype='int32')
adata.var_names_make_unique()
adatas.append(adata)

adata = adatas[0].concatenate(adatas[1],index_unique=None)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Glioblastoma/GSE146026/adata.h5ad"
adata.write(results_file)

####################
#E-MTAB-8107
######################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/E-MTAB-8107')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/E-MTAB-8107')

import re
a_l=[]
for i in files:
	a=re.split(".counts.csv",i)
	a_l.append(a[0])

adatas = []
for x in files:
	adata= sc.read_csv(x)
	adata=adata.T
	adatas.append(adata)

var_names = adatas[0].var_names.union(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.union(adatas[j].var_names)

adata = adatas[0].concatenate(adatas[1:],batch_categories=a_l)
adata.obs['samples']=adata.obs.batch
adata.obs['patients']=adata.obs.batch
adata.obs_names_make_unique()
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/E-MTAB-8107/adata.h5ad"
adata.write(results_file)

###############
#整合
###############

adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/E-MTAB-8107/adata.h5ad')
adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/OV_FTC_GSE154763/adata.h5ad')

adata2=adata2[adata2.obs['tissue']=='T']
adata2.obs['patients']=adata2.obs['patient']
adata2.obs['samples']=adata2.obs['patient']
adata2.obs=adata2.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch','library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster','source', 'tech', 'UMAP1', 'UMAP2'],axis=1)
adata1.obs=adata1.obs.drop(labels=['batch'],axis=1)

#var_names = adata1.var_names.union(adata2.var_names) #

#adata1=re_size_adata(adata1,var_names)
#adata2=re_size_adata(adata2,var_names)

adata_concat = adata1.concatenate(adata2,batch_categories=['E-MTAB-8107', 'GSE154763']) #合并数据集
adata_concat.obs['datasets']=adata_concat.obs['batch']
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)

import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='datasets')  
adata=adata_concat

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'OvarianCancer_integrate_cell.h5ad'
adata_concat.write(results_file)

###############
#细胞类型注释
##############
import scanpy as sc
import celltypist
from celltypist import models
#使用 CellTypist 内置模型分配细胞类型标签
#models.download_models(force_update = True)
#了解模型及其代表的内容。
#models.models_description()
#选择模型
model = models.Model.load(model = 'Pan_Fetal_Human.pkl')

#adata_2000.X.expm1().sum(axis = 1)

sc.pp.normalize_total(adata_concat, target_sum=1e4) ##标准化
sc.pp.log1p(adata_concat)
sc.pp.highly_variable_genes(adata_concat)

predictions1 = celltypist.annotate(adata_concat, model = 'Pan_Fetal_Human.pkl')
adata_concat = predictions1.to_adata()
del predictions1
#############################
#3.5 删除无用细胞
for i in ['EPITHELIUM_II', 'CYCLING_EPITHELIUM','MYOFIBROBLAST','EPITHELIUM_I','HEPATOCYTE_II','HEPATOCYTE_I','ENDOTHELIUM_IV', 'HEPATOCYTE-LIKE', 'VSMC_PERICYTE','MEMP','MESOTHELIUM','VSMC_PERICYTE_III','ENTEROENDOCRINE_I','VSMC_PERICYTE_I','MESENCHYMAL_LYMPHOID_TISSUE_ORGANISER','ENDOTHELIUM_II','DEVELOPING_NEPHRON_I','ENTEROENDOCRINE_II','ENDOTHELIUM_III','VSMC_PERICYTE_II','ENDOTHELIUM_V','LANGERHANS_CELLS','HIGH_MITO','DOUBLETS_FIBRO_ERY','FIBROBLAST_XVII','FIBROBLAST_XI','FIBROBLAST_I','FIBROBLAST_VIII','FIBROBLAST_X','FIBROBLAST_IX','CYCLING_FIBROBLAST_II','FIBROBLAST_XIII','FIBROBLAST_VII','FIBROBLAST_II', 'FIBROBLAST_XII','FIBROBLAST_XV','MOP','ABT(ENTRY)','YS_ERY','MID_ERY','NEURON','CYCLING_YS_ERY','DOUBLET_ENDOTHELIUM_ERYTHROCYTE','LOW_Q_INCONSISTENT','CYCLING_FIBROBLAST_I','FIBROBLAST_XVI','FIBROBLAST_IV','ENDOTHELIUM_I','LATE_ERY','FIBROBLAST_V','KERATINOCYTE','YS_STROMA','SMOOTH_MUSCLE','FIBROBLAST_VI','INTERSTITIAL_CELLS_OF_CAJAL','SKELETAL_MUSCLE','OSTEOCLAST','MUSCLE_SATELLITE','DOUBLET_VSMC_ERYTHROCYTE','FIBROBLAST_III','PLACENTAL_CONTAMINANTS','MELANOCYTE','DEVELOPING_NEPHRON_II','CHONDROCYTE']:
    adata_concat = adata_concat[adata_concat.obs.predicted_labels != i]

adata_concat.obs_names_make_unique()
adata=sc.read_h5ad('OvarianCancer_integrate_cell.h5ad')

adata.obs_names_make_unique()
adata=adata[adata_concat.obs_names]


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'OvarianCancer_Immune_cell.h5ad'
adata.write(results_file)
adata=sc.read_h5ad(results_file)

```


### 

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
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'OvarianCancer_Immune_cell.h5ad'
adata = sc.read_h5ad(results_file)
sc.external.pp.bbknn(adata, batch_key='patients')  

adata.layers["counts"] = adata.X

sc.pp.filter_genes(adata, min_cells=20)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=4000)
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

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata.write(results_file)
```

### annotation

```python

import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

predictions_low = celltypist.annotate(adata_celltypist, model=model_low, majority_voting=True)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "majority_voting"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()

#adata = adata[adata.obs["celltypist_cell_label_fine"] != "Fibroblasts"]

sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

# 使用pivot_table计算每个类别中不同细胞类型的计数
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

# 根据每个类别中占比最大的细胞类型来确定细胞名
celltype_names = pivot_table.idxmax(axis=1)

adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)

adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_ov_celltypist_cell_label_fine')
adata.obs["louvain_1_anno"].value_counts()
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
sc.external.pp.bbknn(adata, batch_key='patients')

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)
sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

# 使用pivot_table计算每个类别中不同细胞类型的计数
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

# 根据每个类别中占比最大的细胞类型来确定细胞名
celltype_names = pivot_table.idxmax(axis=1)

adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)

adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_ov_celltypist_cell_label_fine')

adata.write(results_file)
```



