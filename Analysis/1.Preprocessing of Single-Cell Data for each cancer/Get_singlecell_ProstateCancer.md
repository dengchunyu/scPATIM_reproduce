```python
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
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())

def un_gz(file_name):
    f_name = file_name.replace(".gz", "")
    g_file = gzip.GzipFile(file_name)
    open(f_name, "wb+").write(g_file.read())
    g_file.close() #关闭gzip对象


def un_tar(file_name):
    #untar zip file
    tar = tarfile.open(file_name)
    names = tar.getnames()
    if os.path.isdir(file_name + "_files"):
        pass
    else:
        os.mkdir(file_name + "_files")
    for name in names:
        tar.extract(name, file_name + "_files/")
    tar.close()
#import file as f

('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer')
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer')
for x in files:
	os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/'+ str(x))
	f_n =str(x) + '_RAW.tar'
	un_tar(f_n)
#####################
#GSE237603
#########################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE237603/')
un_tar("GSE237603_RAW.tar")
prefix=['GSM7634626_scaf1511-','GSM7634630_scaf1427-','GSM7634634_scaf1424-','GSM7634627_scaf1426-','GSM7634631_scaf1512-','GSM7634635_scaf1933-','GSM7634628_scaf1510-','GSM7634632_scaf1425-','GSM7634629_scaf1513-','GSM7634633_scaf1385-']
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE237603/GSE237603_RAW.tar_files/'
adatas = []
for x in prefix:
    adata=sc.read_10x_mtx(paths,prefix=x)
    adata.var_names_make_unique()
    strlist = x.split('_')
    adata.obs['patients']= strlist[1]
    adata.obs['samples']= strlist[0]
    adatas.append(adata)
adata = adatas[0].concatenate(adatas[1:],index_unique=None)

adata.write('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE237603/adata.h5ad')


###############
#GSE176031 基因数量太少，细胞数也少
################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031/GSE176031_RAW.tar_files')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031/GSE176031_RAW.tar_files')
for x in files:
	un_gz(x)
adatas = []
ss=[]
for i in ['GSM5353224_PA_PR5186_Pool_1_2_3_S27_L001_dge.txt','GSM5353225_PA_PR5196-1_Pool_1_2_3_S53_L002_dge.txt','GSM5353226_PA_PR5196-2_Pool_1_2_3_S54_L002_dge.txt','GSM5353227_PA_PR5199-193K_Pool_1_2_3_S55_L002_dge.txt','GSM5353232_PA_PR5249_T1_S3_L001_dge.txt','GSM5353236_PA_PR5251_T1_S7_L001_dge.txt','GSM5353237_PA_PR5251_T2_S8_L001_dge.txt']:
	adata= sc.read_csv(i,delimiter='\t')
	adata=adata.T
	strlist = i.split('_')
	ss.append(strlist[0])
	adata.obs['patients']= strlist[2]
	adatas.append(adata)

var_names = adatas[0].var_names.intersection(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.intersection(adatas[j].var_names)

from scipy.sparse import *
import numpy as np

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=ss)

adata.obs['samples']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='samples')
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031/adata.h5ad"
adata.write(results_file)


################
#GSE141445
#################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE141445/GSE141445_RAW.tar_files')
un_gz("GSM4203181_data.raw.matrix.txt.gz")
adata= sc.read_csv('GSM4203181_data.raw.matrix.txt',delimiter='\t')
adata=adata.T

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE141445/adata.h5ad"

adata.obs['patients']= 'GSM4203181'
adata.obs['samples']= 'GSM4203181'

adata.write(results_file)

#################
#GSE137829
################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/GSE137829_RAW.tar_files')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/GSE137829_RAW.tar_files')
for x in files:
	un_gz(x)

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/GSE137829_RAW.tar_files')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/GSE137829_RAW.tar_files')
adatas = []
for x in files:
	X= pd.read_csv(x,delimiter='\t',index_col=1)
	X= X.iloc[:,1:]
	n_obs=len(X.columns)
	n_vars = len(X.index)
	obs = pd.DataFrame(index=X.columns)
	var = pd.DataFrame(index=X.index)
	X=X.T
	adata=ad.AnnData(X,obs=obs,var=var,dtype='float64')
	strlist = x.split('_')
	adata.obs['patients']= strlist[1]
	adata.obs['samples']= strlist[0]
	adata.var_names_make_unique()
	adata.obs_names_make_unique()
	adatas.append(adata)

#var_names = adatas[0].var_names.intersection(adatas[1].var_names)
#for j in range(len(adatas)):
#    len(adatas[j].var_names)
#    var_names = var_names.intersection(adatas[j].var_names)

#for j in range(len(adatas)):
#    adatas[j]=re_size_adata(adatas[j],var_names)

var_names = []
for j in range(len(adatas)):
    var_names.extend(adatas[j].var_names)

List = var_names
List_set = set(List) #List_set是另外一个列表，里面的内容是List里面的无重复 项
from collections import Counter
Clist=Counter(List)
tf_var_names = []
la=round(1/3 * len(adatas))
for item in List_set:
    if Clist[item] > la:
        tf_var_names.append(item)

tf_var_names=set(tf_var_names)
#补足基因的数量
adata_concat=[]
for j in range(len(adatas)):
    adata_concat.append(re_size_adata(adatas[j],tf_var_names))

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=files)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')

#adata.obs['samples']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/adata.h5ad"
adata.write(results_file)
################
#GSE246684
################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE246684/')
X= pd.read_csv('GSE246684_subread_counts_gene_symbol.txt.gz',delimiter='\t',index_col=1)

################
#整合所有数据，并注释细胞

adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE141445/adata.h5ad")
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/adata.h5ad")
adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031/adata.h5ad")
adata4 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE237603/adata.h5ad')

adata = adata1.concatenate(adata2,adata3,adata4,batch_categories=['GSE141445', 'GSE137829','GSE176031','GSE237603'])
adata.obs_names_make_unique()
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.obs['datasets']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'ProstateCancer_integrate_cell.h5ad'
adata.write(results_file)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'ProstateCancer_Immune_cell.h5ad'
adata.write(results_file)
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

sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/')
results_file = 'ProstateCancer_Immune_cell.h5ad'
adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/ProstateCancer_integrate_cell.h5ad")

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
adata.write(results_file)

def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)
adata=reshape_anno(adata)
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_prostate_celltypist_cell_label_fine')
adata.write(results_file)

adata = adata[adata.obs["louvain_1_anno"] != "Fibroblasts"]
adata = adata[adata.obs["louvain_1_anno"] != "Endothelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_celltypist_cell_label_fine')

adata.write(results_file)

adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/ProstateCancer_Immune_cell.h5ad")
adata = adata[adata.obs["celltypist_cell_label_fine"] != "Epithelial cells"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata=reshape_anno(adata)

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_celltypist_cell_label_fine')
adata = sc.read_h5ad(results_file)
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Fibroblasts"]
adata = adata[adata.obs["louvain_1_anno"] != "Endothelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Double-positive thymocytes"]

adata.obs["louvain_1_anno"].value_counts()

sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)
adata=reshape_anno(adata)

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_prostate_celltypist_cell_label_fine')
adata.write(results_file)


```