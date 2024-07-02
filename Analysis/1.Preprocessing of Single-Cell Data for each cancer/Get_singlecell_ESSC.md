```python
#ESSC
'/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer'
import tarfile
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import gzip
import numpy as np
import warnings
import gc
import anndata as ad
from anndata import AnnData
from scipy.sparse import *
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())
import anndata as ad
###1）解压
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

#########################
#GSE145370
#############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370')
un_tar("GSE145370_RAW.tar")
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files')

for i in files:
	un_gz(i)

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files/matrix')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files/matrix')
for i in files:
	un_tar(i)

filesn = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files/matrix')
paths1='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/GSE145370_RAW.tar_files/matrix'

#2）循环读取mtx文件
adatas = []
for preFix in filesn:
    paths2=paths1+'/'+preFix+'/filtered_feature_bc_matrix'
    adata = sc.read_10x_mtx(paths2)
    strlist = preFix.split('_')
    adata.obs['samples']= strlist[0]
    adata.obs['patients']= strlist[1]
    adata.var_names_make_unique()
    adatas.append(adata)

#3）整合以上样本
adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=filesn)
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
adata.var=adata.var.drop(labels=['feature_types'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/adata.h5ad"
adata.write(results_file)

####################
#GSE160269
###################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269')

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269')

for i in files:
	un_gz(i)

adatas = []
for i in ['Tcell','Myeloid','Bcell']:
    mtx_file_counts='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269/GSE160269_UMI_matrix_' + i +'.txt'
    adata = sc.read_csv(mtx_file_counts,delimiter=' ')
    adata = adata.T
    adata.var_names_make_unique()
    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=['Tcell','Myeloid','Bcell'])

ss=[]
for x in adata.obs_names:
    ss.append(x.split('-')[0])

new_li1 = list(set(ss))
adata.obs['samples']= ss
adata.obs['patients']= ss
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

adata=adata[adata.obs.samples!='P127N']
adata=adata[adata.obs.samples!='P128N']
adata=adata[adata.obs.samples!='P126N']

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269/adata.h5ad"
adata.write(results_file)

##############
#GSE156728
################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE156728')
adata=sc.read_h5ad('adata.h5ad')
adata.obs['patients']=adata.obs['patient']
adata.obs['samples']=adata.obs['patient']
adata.obs=adata.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE156728/adata.h5ad"
adata.write(results_file)
##############
#GSE154763 删除，非count数据
################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE154763')
adata=sc.read_h5ad('adata.h5ad')
adata=adata[adata.obs['tissue']=='T']
adata.obs['patients']=adata.obs['patient']
adata.obs['samples']=adata.obs['patient']
adata.obs=adata.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch', 'library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster', 'source', 'tech', 'UMAP1', 'UMAP2'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE154763/adata.h5ad"
adata.write(results_file)

##################
#合并数据
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
import scipy
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE156728/adata.h5ad")
adata1.obs['patients']=adata1.obs['patient']
adata1.obs['samples']=adata1.obs['patient']
adata1.obs=adata1.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/adata.h5ad")
#adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE154763/adata.h5ad")
adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269/adata.h5ad")

var_names = adata1.var_names.intersection(adata2.var_names)
#var_names = var_names.union(adata3.var_names)
var_names = var_names.intersection(adata4.var_names)

adata = adata1.concatenate(adata2,adata4, batch_categories=['GSE156728','GSE145370','GSE160269']) #合并数据集
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.obs['datasets']=adata.obs.batch
sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='ESCC'
n_comps=50
sc.pl.pca(adata, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')

adata.layers["counts"] = adata.X
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Esophagealcancer_integrate_cell.h5ad'
adata.write(results_file)
adata=sc.read_h5ad(results_file)

import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")
#model_high = models.Model.load(model="Immune_All_High.pkl")
#predictions_high = celltypist.annotate(adata, model=model_high, majority_voting=True)
#predictions_high_adata = predictions_high.to_adata()
adata.obs_names_make_unique()

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

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_liver_celltypist_cell_label_fine')
adata.obs["louvain_1_anno"].value_counts()

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
results_file = 'EsophagealCancer_Immune_cell.h5ad'
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


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='EsophagealCancer_Immune_cell.h5ad'

adata = sc.read_h5ad('EsophagealCancer_Immune_cell.h5ad')

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

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)

adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

adata.layers["counts"] = adata.X
sc.pp.filter_genes(adata, min_cells=20)
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

#adata = adata[adata.obs["celltypist_cell_label_fine"] != "Fibroblasts"]

sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

# 使用pivot_table计算每个类别中不同细胞类型的计数
pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

# 根据每个类别中占比最大的细胞类型来确定细胞名
celltype_names = pivot_table.idxmax(axis=1)

adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)

adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_escc_celltypist_cell_label_fine')
adata.obs["louvain_1_anno"].value_counts()
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_escc_celltypist_cell_label_fine')

adata.write(results_file)
```