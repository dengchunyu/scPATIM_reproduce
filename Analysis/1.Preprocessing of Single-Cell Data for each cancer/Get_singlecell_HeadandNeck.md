```python
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck
from scipy.sparse import *
import anndata as ad

'/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck'
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck')
for i in files:
	f='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/'+i+'/'+i+'_RAW.tar'
	un_tar(f)

############
#GSE139324
###########
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/GSE139324_RAW.tar_files')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/GSE139324_RAW.tar_files')
for i in files:
	un_gz(i)

a=np.arange(27)
a=a[1:27]
b=a*2+9

paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/GSE139324_RAW.tar_files'
adatas = []
pn=[]
for i in range(26):
    prefix='GSM41381'+str(b[i])+'_HNSCC_'+str(a[i])+'_TIL_'
    p='GSM41381'+str(b[i])
    pn.append(p)
    adata = sc.read_10x_mtx(paths,prefix=prefix)
    adata.var_names_make_unique()
    adatas.append(adata)



adata = adatas[0].concatenate(adatas[1:],batch_categories=pn)
adata.obs['samples']=adata.obs.batch
adata.obs['patients']=adata.obs.batch
adata.obs['subCancerType']='HNSCC'

adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/adata.h5ad"
adata.write(results_file)

###########
#GSE164690
##########

files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/GSE164690_RAW.tar_files')

sample_ids=['GSM5017022','GSM5017025','GSM5017027','GSM5017029','GSM5017031','GSM5017034','GSM5017037','GSM5017040','GSM5017043','GSM5017046','GSM5017049','GSM5017052','GSM5017055','GSM5017058','GSM5017061','GSM5017064','GSM5017067','GSM5017070']
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/GSE164690_RAW.tar_files')


current_directory = os.getcwd()  # 获取当前工作目录

# 使用列表推导式列出以指定样本ID为开头的文件
matching_files = [filename for filename in os.listdir(current_directory) if any(filename.startswith(sample_id) for sample_id in sample_ids)]



import re
a_l=[]
for i in files:
    bc=re.compile("barcodes.tsv.gz")
    bci=bc.search(i)
    if bci!=None:
        a=re.split("barcodes.tsv.gz",i)
        a_l.append(a[0])

def de_duplication(list):
    dedup_list = []
    for word in list:
        if not word in dedup_list:
            dedup_list.append(word)
    return dedup_list

prefile=de_duplication(a_l)
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/GSE164690_RAW.tar_files'
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/GSE164690_RAW.tar_files')

adatas = []
for preFix in prefile:
    un_gz(preFix+'matrix.mtx.gz')
    os.remove(preFix+'matrix.mtx.gz')
    un_gz(preFix+'barcodes.tsv.gz')
    os.remove(preFix+'barcodes.tsv.gz')
    un_gz(preFix+'features.tsv.gz')
    os.remove(preFix+'features.tsv.gz')
    os.rename(preFix+'features.tsv', preFix+'genes.tsv')
    adata = sc.read_10x_mtx(paths,prefix=preFix)
    adata.var_names_make_unique()
    adatas.append(adata)



for j in range(len(adatas)):
    len(adatas[j].var_names)

adata = adatas[0].concatenate(adatas[1:],batch_categories=prefile[0:len(adatas)])

ss=[]
pp=[]
for i in adata.obs.batch:
    a=re.split("_",i)
    ss.append(a[0])
    pp.append(a[1])

adata.obs['samples']=ss
adata.obs['patients']=pp
adata.obs['subCancerType']='PBL'
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/adata.h5ad"
adata.write(results_file)

#############
#
#adata0 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/adata.h5ad")
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/adata.h5ad")
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/adata.h5ad")

#var_names = adata1.var_names.union(adata2.var_names) #
#adata1=re_size_adata(adata1,var_names)
#adata2=re_size_adata(adata2,var_names)

adata = adata1.concatenate(adata2,batch_categories=['GSE139324', 'GSE164690'])
#adata = adata.concatenate(adata3,index_unique=None)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.obs['datasets']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'HeadandNeck_integrate_cell.h5ad'
adata.write(results_file)


import scanpy as sc
import celltypist
from celltypist import models
#使用 CellTypist 内置模型分配细胞类型标签
models.download_models(force_update = True)
#了解模型及其代表的内容。
models.models_description()
#选择模型
model = models.Model.load(model = 'Pan_Fetal_Human.pkl')

#adata_2000.X.expm1().sum(axis = 1)

sc.pp.normalize_total(adata, target_sum=1e4) ##标准化
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
predictions1 = celltypist.annotate(adata, model = 'Pan_Fetal_Human.pkl')
adata = predictions1.to_adata()
adata.obs.predicted_labels.unique()[1:10]
del predictions1
for i in ['EPITHELIUM_II', 'CYCLING_EPITHELIUM','MYOFIBROBLAST','EPITHELIUM_I','HEPATOCYTE_II','HEPATOCYTE_I','ENDOTHELIUM_IV', 'HEPATOCYTE-LIKE', 'VSMC_PERICYTE','MEMP','MESOTHELIUM','VSMC_PERICYTE_III','ENTEROENDOCRINE_I','VSMC_PERICYTE_I','MESENCHYMAL_LYMPHOID_TISSUE_ORGANISER','ENDOTHELIUM_II','DEVELOPING_NEPHRON_I','ENTEROENDOCRINE_II','ENDOTHELIUM_III','VSMC_PERICYTE_II','ENDOTHELIUM_V','LANGERHANS_CELLS','HIGH_MITO','DOUBLETS_FIBRO_ERY','FIBROBLAST_XVII','FIBROBLAST_XI','FIBROBLAST_I','FIBROBLAST_VIII','FIBROBLAST_X','FIBROBLAST_IX','CYCLING_FIBROBLAST_II','FIBROBLAST_XIII','FIBROBLAST_VII','FIBROBLAST_II', 'FIBROBLAST_XII','FIBROBLAST_XV','MOP','ABT(ENTRY)','YS_ERY','MID_ERY','NEURON','CYCLING_YS_ERY','DOUBLET_ENDOTHELIUM_ERYTHROCYTE','LOW_Q_INCONSISTENT','CYCLING_FIBROBLAST_I','FIBROBLAST_XVI','FIBROBLAST_IV','ENDOTHELIUM_I','LATE_ERY','FIBROBLAST_V','KERATINOCYTE','YS_STROMA','SMOOTH_MUSCLE','FIBROBLAST_VI','INTERSTITIAL_CELLS_OF_CAJAL','SKELETAL_MUSCLE','OSTEOCLAST','MUSCLE_SATELLITE','DOUBLET_VSMC_ERYTHROCYTE','FIBROBLAST_III','PLACENTAL_CONTAMINANTS','MELANOCYTE','DEVELOPING_NEPHRON_II','CHONDROCYTE','MYOFIBROBLAST_I']:
    adata = adata[adata.obs.predicted_labels != i]


adata_concat=sc.read_h5ad(results_file)
adata_concat=adata_concat[adata.obs_names]

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')

results_file = 'HeadandNeck_Immune_cell.h5ad'
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


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')

sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='HeadandNeck_Immune_cell.h5ad'

adata = sc.read_h5ad(results_file)

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

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata.write(results_file)
```

### annotation

```python

import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")

duplicated_obs_names = adata.obs_names.duplicated(keep='first')
adata = adata[~duplicated_obs_names]

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

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_hnsc_celltypist_cell_label_fine')
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
adata.write(results_file)
```