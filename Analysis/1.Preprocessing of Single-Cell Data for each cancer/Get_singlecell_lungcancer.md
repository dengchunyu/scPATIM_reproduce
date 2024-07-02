## lung cancer 数据的读取

```python
import tarfile
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import numpy as np
import warnings
from anndata import AnnData
import matplotlib.pyplot as pt
import gc
import gzip
import scipy
from scipy.sparse import *
    
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
```


### GSE131907

```python
anno_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"

matrix_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds.gz"

un_gz(anno_file)
un_gz(matrix_file)

anno_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt"

matrix_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds"

print('Reading metadata')
anno_data = pd.read_table(anno_file,index_col=0)
anno_data['Sample_Origin'].unique()
anno_data['Cell_type'].unique()

print('Reading count matrix')
matrix_data = sc.read_text(matrix_file)

print('Transposing matrices')
matrix_data = matrix_data.T
# Append metadata
print('Appending metadata')
matrix_data.obs = anno_data.loc[matrix_data.obs_names]
#
out_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata.h5ad"

matrix_data.write(out_file)
adata1 = sc.read_h5ad(out_file)
adata1 = adata1[adata1.obs.Sample_Origin != 'nLung']
adata1 = adata1[adata1.obs.Sample_Origin != 'nLN']
adata1 = adata1[adata1.obs.Sample_Origin != 'mBrain']
adata1 = adata1[adata1.obs.Sample_Origin != 'PE']
adata1.write(out_file)
adata1.obs['patients']=adata1.obs['Sample']
adata1.obs['samples']=adata1.obs['Sample']
adata1.write("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata1.h5ad")
adata1.obs=adata1.obs.drop(labels=['Barcode', 'Sample', 'Sample_Origin','Cell_type.refined', 'Cell_subtype'],axis=1)
adata1.write("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata.h5ad")
```


### Leader_et_al

```python
cell_metadata = pd.read_csv('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al/cell_metadata.csv')
cell_metadata.cell_ID.to_csv( 'barcodes.tsv', sep='\t',index=False,header=False)

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al')
adata2 = sc.read_10x_mtx("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al/",cache=True)
gc.collect()
metadata = pd.read_csv('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al/cell_metadata.csv', index_col=0)
adata2.obs = metadata.loc[adata2.obs_names]

##删除正常样本
for i in [37,39,41,44,46,49,50,54,55,56,60,61,63,96,98,100,104,108,112,115,310,312,314,316,319,323,325,327,331,333,335,337,339,341,343,399,413,449,453,457,461,465,469,473,480,64,74,75,82,86,90,91,360,361,362,363,366,370,376,379]:
	adata2 = adata2[adata2.obs.sample_ID != i]

adata2.obs['patients']=adata2.obs['sample_ID']
adata2.obs['samples']=adata2.obs['sample_ID']
adata2.obs['datasets']='Leader_et_al'
adata2.obs=adata2.obs.drop(labels=['sample_ID', 'cluster_ID'],axis=1)

out_file = 'adata_all.h5ad'
adata2.write(out_file)

##删除肿瘤细胞
for i in [3,4,6,12,15,21,22,24,26,27,60]:
	adata2 = adata2[adata2.obs.cluster_ID != i]
        
out_file = 'adata.h5ad'
adata2.write(out_file)

```

### E-MTAB-6653 fastq从头处理数据

```python
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653')

filesn = ["BT1375_S5_L004","BT1376_S6_L005","BT1377_S7_L006","scrBT1426_hg19_S12_L006","scrBT1426_hg19_S12_L006","scrBT1427_hg19_S13_L007","scrBT1430m_S0_L002","scrBT1431m_S0_L003","scrBT1432m_S0_L004"]
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653'
#2）循环读取mtx文件
adatas = []
for preFix in filesn:
    path= paths+"/"+ preFix+ "/outs/filtered_feature_bc_matrix"
    adata = sc.read_10x_mtx(path)
    strlist = preFix.split('_')
    adata.obs['samples']= strlist[0]
    adata.obs['patients']= strlist[0]
    adata.var_names_make_unique()
    adatas.append(adata)

var_names = adatas[0].var_names.intersection(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.intersection(adatas[j].var_names)

#3）整合以上样本
adata = adatas[0].concatenate(adatas[1:],index_unique=None)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='samples')

adata.obs=adata.obs.drop(labels=['batch'],axis=1)
adata.var=adata.var.drop(labels=['feature_types'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653/adata.h5ad"
adata.write(results_file)
```

### 从上面的数据中分别输出数据矩阵方便R读取

```python
#合并
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')

adata1 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata.h5ad')

adata2 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al/adata_all.h5ad')

#adata2是我们需要计算的
adata3 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653/adata.h5ad')

adata2.var_names=adata2.var.gene_ids
adata2.obs.patients
ps=[]
for x in adata2.obs.patients:
    x="Patient_"+str(x)
    ps.append(x)

adata2.obs.patients= ps
adata2.obs.samples= ps

adata_concat = adata1.concatenate(adata2,adata3, batch_categories=['GSE131907', 'Leader_et_al','E-MTAB-6653'])
adata_concat.obs_names_make_unique()
adata_concat.obs['datasets']=adata_concat.obs['batch']

import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')  
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Lung_integrate_cell.h5ad'
adata_concat.write(results_file)

adata2.obs_names_make_unique()
adata1.obs_names_make_unique()
adata3.obs_names_make_unique()

var_names = adata1.var_names.intersection(adata2.var_names)
var_names = var_names.intersection(adata3.var_names)
adata1=adata1[:,var_names]
adata2=adata2[:,var_names]
adata3=adata3[:,var_names]

sc.pp.pca(adata1)
sc.pp.neighbors(adata1)
sc.tl.umap(adata1)
#sc.pl.umap(adata4, color='Coarse_Cell_Annotations',legend_loc='on data',legend_fontsize='xx-small',save='pancreas_SCP1644_umap.pdf')

sc.tl.ingest(adata2, adata1, obs='Cell_type')
sc.tl.ingest(adata3, adata1, obs='Cell_type')

#adata1 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata.h5ad')

for i in ['Epithelial cells','Endothelial cells','Fibroblasts']:
    adata2 = adata2[adata2.obs.Cell_type != i]
    adata1 = adata1[adata1.obs.Cell_type != i]
    adata3 = adata3[adata3.obs.Cell_type != i]

adata1.obs_names_make_unique()
adata2.obs_names_make_unique()
adata3.obs_names_make_unique()

adata_concat = adata1.concatenate(adata2,adata3, batch_categories=['GSE131907', 'Leader_et_al','E-MTAB-6653']) #合并数据集
adata_concat.obs=adata_concat.obs.drop(labels=['Cell_type'],axis=1)

import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')  
adata_concat.obs['datasets']=adata_concat.obs['batch']
#adata_concat.obs=adata_concat.obs.drop(labels=['patients', 'samples'],axis=1)

adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)

results_file = 'Lung_Immune_cell.h5ad'
adata_concat.write(results_file)

color_map='CMRmap'
cancer='Lung'
n_comps=50
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')


##################
#3.2 过滤基因和细胞
#sc.pp.filter_cells(adata_concat, min_genes=200)
#sc.pp.filter_genes(adata_concat, min_cells=3)
#mito_genes=adata_concat.var_names.str.startswith('MT-')
#adata_concat.obs['percent_mito']=np.sum(adata_concat[:,mito_genes].X,axis=1).A1/np.sum(adata_concat.X,axis=1).A1
#adata_concat.obs['n_counts']=adata_concat.X.sum(axis=1).A1
##过滤线粒体基因比例 > 5% 和基因总数 >2500 的细胞。
#adata_concat = adata_concat[adata_concat.obs.n_genes < 7000, :]
#adata_concat = adata_concat[adata_concat.obs.percent_mito < 0.5, :]
gc.collect()
#########################
#3.3 预处理
#sc.pp.normalize_total(adata_concat, target_sum=1e4) ##标准化
#sc.pp.log1p(adata_concat)
#sc.pp.highly_variable_genes(adata_concat)
#results_file = 'Lung_Immune_cell.h5ad'
#adata_concat.write(results_file)
```

## 单独处理

### Quality control
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

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
results_file = 'LungCancer_Immune_cell.h5ad'
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')
adata = sc.read_h5ad('Lung_integrate_cell.h5ad')
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

adata.obs.mt_outlier.value_counts()
print(f"Total number of cells: {adata.n_obs}")
#Total number of cells: 193699
adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
#1156224

adata.layers["counts"] = adata.X

sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")
adata.write(results_file)

adata = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/LungCancer_Immune_cell.h5ad')

```

### normalization

```python
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)
adata.obs.patients= "P"+ str(adata.obs.samples.value())
```

### annotation

```python
def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)
model_low = models.Model.load(model="Immune_All_Low.pkl")

adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform


predictions_low = celltypist.annotate(adata_celltypist, model=model_low, majority_voting=True)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "majority_voting"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()
adata=reshape_anno(adata)
for x in ['Endothelial cells','Epithelial cells','Fibroblasts']:
    adata = adata[adata.obs["louvain_1_anno"] != x]

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata=reshape_anno(adata)
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_lung_celltypist_cell_label_fine')

adata.write(results_file)
```