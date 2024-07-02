

```python
import tarfile
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import infercnvpy as cnv
import numpy as np
import warnings
import gc
import gzip
from scipy.sparse import *
import anndata as ad
from anndata import AnnData
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
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

##############
#GSE123139
###############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139')
un_tar('GSE123139_RAW.tar')

files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/GSE123139_RAW.tar_files')

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/GSE123139_RAW.tar_files')
for x in files:
	un_gz(x)

files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/GSE123139_RAW.tar_files')
adatas = []
for x in files:
    adata=sc.read_csv(x,delimiter='\t')
    adata=adata.T
    strlist = x.split('.')
    strlist2=strlist[0].split('_')
    adata.obs['patients']= strlist2[0]
    adata.obs['samples']= strlist2[1]
    adatas.append(adata)

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=files)
#adata.obs['samples']=adata.obs.batch
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')

adata.obs=adata.obs.drop(labels=['batch'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/adata.h5ad"
adata.write(results_file)


############
#GSE115978
##############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE115978')
un_gz('GSE115978_counts.csv.gz')
un_gz('GSE115978_cell.annotations.csv.gz')

adata=sc.read_csv('GSE115978_counts.csv')
adata=adata.T
meta=pd.read_csv('GSE115978_cell.annotations.csv')
adata.obs=meta
adata.obs=adata.obs.drop(labels=['cells', 'treatment.group', 'Cohort', 'no.of.genes', 'no.of.reads'],axis=1)
adata.obs['patients']=adata.obs.samples

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE115978/adata.h5ad"
adata.write(results_file)


###################
#整合

adata1=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/adata.h5ad")
adata1.var_names_make_unique()
adata2=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE115978/adata.h5ad")
adata2.var_names_make_unique()
var_names = adata1.var_names.intersection(adata2.var_names) #
adata1=adata1[:,var_names]
adata2=adata2[:,var_names]
sc.pp.pca(adata2)
sc.pp.neighbors(adata2)
sc.tl.umap(adata2)
sc.tl.ingest(adata1, adata2, obs='cell.types')

adata_concat = adata1.concatenate(adata2, batch_categories=['GSE123139', 'GSE115978'])
adata_concat.obs['datasets']=adata_concat.obs.batch
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)

import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='datasets')  

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Melanoma_integrate_cell.h5ad'
adata_concat.write(results_file)

adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'Mal']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'Endo.']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'CAF']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != '?']
adata_concat.obs=adata_concat.obs.drop(labels=['cell.types'],axis=1)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Melanoma_Immune_cell.h5ad'
from scipy.sparse import *
import numpy as np
adata_concat.X = csc_matrix(adata_concat.X,dtype='float64')
adata_concat.write(results_file)





adata1_tmp=adata1

adata2=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE115978/adata.h5ad")
adata1=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/adata.h5ad")

var_names = adata1.var_names.union(adata2.var_names) #
adata1=re_size_adata(adata1,var_names)
adata2=re_size_adata(adata2,var_names)
adata1.obs['cell.types']=adata1_tmp.obs['cell.types']
#['Mal', 'Macrophage', '?', 'Endo.', 'T.CD4', 'CAF', 'T.CD8',
 #      'T.cell', 'NK', 'B.cell']
del adata1_tmp
adata_concat = adata1.concatenate(adata2, batch_categories=['GSE123139', 'GSE115978'])
adata_concat.obs['datasets']=adata_concat.obs.batch
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)

import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='datasets')  

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Melanoma_integrate_cell.h5ad'
adata_concat.write(results_file)

adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'Mal']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'Endo.']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != 'CAF']
adata_concat = adata_concat[adata_concat.obs['cell.types'] != '?']
adata_concat.obs=adata_concat.obs.drop(labels=['cell.types'],axis=1)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Melanoma_Immune_cell.h5ad'
adata_concat.write(results_file)
```


## 单独处理

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

sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/')
results_file = 'Melanoma_Immune_cell.h5ad'

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
adata.layers["counts"] = adata.X

sc.pp.filter_genes(adata, min_cells=20)
print(f"Number of genes after cell filter: {adata.n_vars}")
adata.write(results_file)
adata = sc.read_h5ad(results_file)

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

model_low = models.Model.load(model="Immune_All_Low.pkl")
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

predictions_low = celltypist.annotate(adata_celltypist, model=model_low, majority_voting=True)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "majority_voting"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()

sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

celltype_names = pivot_table.idxmax(axis=1)

adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)

adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_melanoma_celltypist_cell_label_fine')
adata.obs["louvain_1_anno"].value_counts()
adata = adata[adata.obs["louvain_1_anno"] != "Double-positive thymocytes"]
adata.write(results_file)
```

```python

sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

```