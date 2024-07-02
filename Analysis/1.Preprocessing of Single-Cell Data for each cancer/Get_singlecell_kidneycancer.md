```python
from scipy.sparse import *
import anndata as ad
import gzip
import tarfile

#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer
def un_tar(file_name):
    tar = tarfile.open(file_name)
    names = tar.getnames()
    if os.path.isdir(file_name + "_files"):
        pass
    else:
        os.mkdir(file_name + "_files")
    for name in names:
        tar.extract(name, file_name + "_files/")
    tar.close()
####################
#GSE222703
####################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703')
un_tar("GSE222703_RAW.tar")
prefix=['GSM6929206_p022_Tumoral_','GSM6929208_p027_Tumoral_','GSM6929210_p029_Tumoral_']
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703/GSE222703_RAW.tar_files')

paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703/GSE222703_RAW.tar_files/'
adatas = []
for x in prefix:
    adata=sc.read_10x_mtx(paths,prefix=x)
    adata.var_names_make_unique()
    strlist = x.split('_')
    adata.obs['patients']= strlist[1]
    adata.obs['samples']= strlist[0]
    adatas.append(adata)
adata = adatas[0].concatenate(adatas[1:],index_unique=None)
sc.pp.filter_cells(adata, min_counts=100)
#sc.pp.filter_cells(adata, max_counts=5000)
adata.write('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703/adata.h5ad')
##################
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/
###################

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813')
un_tar("GSE202813_RAW.tar")
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/GSE202813_RAW.tar_files/')

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/GSE202813_RAW.tar_files/')
adatas = []
for i in files:
    mtx_file_counts="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/GSE202813_RAW.tar_files/" + i
    adata = sc.read_csv(mtx_file_counts)
    adata = adata.T
    strlist = i.split('-')
    adata.obs['samples']= strlist[0]
    adata.obs['patients']= strlist[1]
    adata.var_names_make_unique()
    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:],index_unique=None)
adata.write('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/adata.h5ad')

################
#GSE121638
#################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638')
#1)解压缩
un_tar("GSE121638_RAW.tar")
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/GSE121638_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/GSE121638_RAW.tar_files')
for i in files:
	un_gz(i)
prefix=['GSM3440846_GU0744_T_','GSM3440845_GU0715_T_','GSM3440844_GU0700_T_']

adatas = []
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/GSE121638_RAW.tar_files'
ss=[]
for x in prefix:
    adata=sc.read_10x_mtx(paths,prefix=x)
    adata.var_names_make_unique()
    strlist = x.split('_')
    ss.append(strlist[0])
    adata.obs['patients']= strlist[1]
    adatas.append(adata)

var_names = adatas[0].var_names.union(adatas[1].var_names)
var_names = var_names.union(adatas[2].var_names)


adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=ss)
adata.obs['samples']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/adata.h5ad"
adata.write(results_file)



adata_concat = adata1.concatenate(adata3,batch_categories=["GSE156728","GSE121638"]) #合并数据集
adata_concat.obs['datasets']=adata_concat.obs.batch
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='datasets')

results_file = 'KidneyCancer_Immune_cell.h5ad'
adata_concat.write(results_file)
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
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/')

results_file = 'KidneyCancer_Immune_cell.h5ad'
adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/RC_GSE156728/adata.h5ad')

adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/KIDNEY_GSE154763/adata.h5ad')
adata3=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/adata.h5ad')
adata4=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703/adata.h5ad')
adata5=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/adata.h5ad')

adata= adata1.concatenate(adata2,adata3,adata4,adata5,batch_categories=["GSE156728","GSE154763","GSE121638","GSE222703","GSE202813"]) #合并数据集
adata.obs['datasets']=adata.obs.batch

sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='datasets')
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

#results_file=os.path.join(BaseDirectory, 'adata_pp_immune.h5ad')
adata.write('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/KidneyCancer_integrate_cell.h5ad')
adata=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/KidneyCancer_integrate_cell.h5ad')
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

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

adata.layers["counts"] = adata.X

sc.pp.filter_genes(adata, min_cells=20)
adata.write(results_file)

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
adata=sc.read_h5ad(results_file)
adata.obs_names_make_unique()
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

def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)

adata=reshape_anno(adata)
adata.obs["louvain_1_anno"].value_counts()

#adata = adata[adata.obs["louvain_1_anno"] != "Fibroblasts"]
#adata = adata[adata.obs["louvain_1_anno"] != "Endothelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_kidney_celltypist_cell_label_fine')

adata.write(results_file)
adata=sc.read_h5ad(results_file)
```