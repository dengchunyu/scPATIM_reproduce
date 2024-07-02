### GSE193581
```python
import tarfile
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import gzip
import infercnvpy as cnv
import numpy as np
import warnings
import gc
import anndata as ad
from anndata import AnnData
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
def un_gz(file_name):
    f_name = file_name.replace(".gz", "")
    g_file = gzip.GzipFile(file_name)
    open(f_name, "wb+").write(g_file.read())
    g_file.close() #关闭gzip对象


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

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/')
#1)解压缩
un_tar("GSE193581_RAW.tar")


files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_RAW.tar_files/')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_RAW.tar_files/')
for i in files:
	un_gz(i)
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_RAW.tar_files/')

adatas = []
for i in files:
    mtx_file_counts="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_RAW.tar_files/" + i
    adata = sc.read_text(mtx_file_counts,delimiter='\t')
    adata = adata.T
    strlist = i.split('_')
    adata.obs['samples']= strlist[0]
    adata.obs['patients']= strlist[1]
    adata.var_names_make_unique()
    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:],index_unique=None)
metadata = pd.read_csv('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_celltype_annotation.txt.gz',sep='\t')

adata.obs['id']=adata.obs_names
metadata['id']=metadata.index
adata.obs = pd.merge(adata.obs, metadata, on='id', how='left')


os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/')
adata = adata[adata.obs['celltype'].notna()]
results_file = 'GSE193581.h5ad'
adata.write(results_file)

adata=sc.read_h5ad(results_file)

#adata.obs['celltype'] = adata.obs['celltype'].astype('category')
selected_cell_types = ['T cell', 'Myeloid cell', 'NK cell', 'B cell']
adata = adata[adata.obs['celltype'].isin(selected_cell_types)]
results_file = 'GSE193581_Immune_cell.h5ad'
#adata.write(results_file)
#adata=sc.read_h5ad(results_file)

##删除样本
adata = adata[adata.obs["patients"] != "NORM03"]
adata = adata[adata.obs["patients"] != "NORM18"]
adata = adata[adata.obs["patients"] != "NORM19"]
adata = adata[adata.obs["patients"] != "NORM20"]
adata = adata[adata.obs["patients"] != "NORM07"]
adata = adata[adata.obs["patients"] != "NORM21"]

sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.write(results_file)


adata.layers["counts"] = adata.X

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

adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

sc.pp.filter_genes(adata, min_cells=20)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)

import seaborn as sns
import matplotlib.pyplot as plt

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)
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

sc.pl.umap(adata,color=["samples","louvain_res1","louvain_1_anno"],save='_celltypist_cell_label_fine')

adata.write(results_file)
```

## GSE154763合并GSE156728
```python
###############

adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/Thyroid_GSE154763/adata.h5ad')
adata1=adata1[adata1.obs['tissue']=='T']
adata1.obs['patients']=adata1.obs['patient']
adata1.obs['samples']=adata1.obs['patient']
adata1.obs=adata1.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch', 'library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster', 'source', 'tech', 'UMAP1', 'UMAP2'],axis=1)

adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/THCA_GSE156728/adata.h5ad')

adata2.obs['patients']=adata2.obs['patient']
adata2.obs['samples']=adata2.obs['patient']
adata2.obs=adata2.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

adata3=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/GSE193581/GSE193581_Immune_cell.h5ad')

adata_concat = adata2.concatenate(adata3,batch_categories=["GSE156728","GSE193581"]) #合并数据集
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')


adata_concat.obs['datasets']=adata_concat.obs.batch
color_map='CMRmap'
cancer='Thyroid'
n_comps=50
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/')
sc.pl.pca(adata_concat, color='samples', color_map=color_map, save='_' + str(cancer))

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/')
results_file = 'ThyroidCancer_immune_cell.h5ad'
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


BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/'
os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='ThyroidCancer_Immune_cell.h5ad'
adata=sc.read_h5ad(results_file)
adata.layers["counts"] = adata.X

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

adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

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

import seaborn as sns
import matplotlib.pyplot as plt

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata.write(results_file)
```

### annotation

```python
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ThyroidCancer/'
os.chdir(BaseDirectory)
adata=sc.read_h5ad(results_file)
import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")


predictions_low = celltypist.annotate(adata, model=model_low, majority_voting=True)
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

sc.pl.umap(adata,color=["samples","louvain_res1","louvain_1_anno"],save='_thy_celltypist_cell_label_fine')
sc.pl.umap(adata,color=["samples","patients"],save='_thy_patients')

for x in ['ATC08','ATC09','ATC10','ATC11','ATC12','ATC13','ATC14','ATC15','ATC16','ATC17']:
    adata = adata[adata.obs["patients"] != x]

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)
sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
celltype_names = pivot_table.idxmax(axis=1)
adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_thy_celltypist_cell_label_fine')

adata.write(results_file)


```