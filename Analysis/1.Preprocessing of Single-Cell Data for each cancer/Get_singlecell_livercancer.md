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

from anndata import AnnData
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())
# load adata
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer')

for i in ['GSE138709']:
	f='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/'+i+'/'+i+'_RAW.tar'
	un_tar(f)
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
    #由于解压后是许多文件，预先建立同名文件夹
    for name in names:
        tar.extract(name, file_name + "_files/")
    tar.close()
import scipy
from scipy.sparse import *
import numpy as np

#import file as f
######################
#GSE151530
###################
'/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530'

files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530')

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530')
for x in files:
	un_gz(x)
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530'
adata = sc.read_10x_mtx(paths,prefix='GSE151530_')
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530/adata.h5ad"

metadata = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530/GSE151530_Info.txt", sep="\t",index_col=0)
adata.obs = metadata#.loc[adata.obs_names]

adata.obs['samples']=adata.obs.Sample
adata.obs['patients']=adata.obs_names
adata.obs=adata.obs.drop(labels=['Sample', 'Cell'],axis=1)
adata.write(results_file)

################
#GSE140228
###############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228/')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228')

for x in files:
	un_gz(x)

paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228'
adata = sc.read_mtx('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228/GSE140228_UMI_counts_Droplet.mtx')
adata=adata.T
genes= pd.read_csv('GSE140228_UMI_counts_Droplet_genes.tsv',delimiter='\t')
adata.var=genes
adata.var.index=genes.SYMBOL
adata.var_names=genes.SYMBOL

meta= pd.read_csv('GSE140228_UMI_counts_Droplet_cellinfo.tsv',delimiter='\t')
barcodes=pd.read_csv('GSE140228_UMI_counts_Droplet_barcodes.tsv',header=None)
adata.obs_names=barcodes
adata.obs=meta
adata.obs['samples']=adata.obs.Sample
adata.obs['patients']=adata.obs.Donor
adata.obs.Tissue_sub.unique()
adata.obs['Tissue_sub'] = adata.obs['Tissue_sub'].astype('category')

adata = adata[adata.obs['Tissue_sub'] != 'Ascites']
adata = adata[adata.obs['Tissue_sub'] != 'Blood']
adata = adata[adata.obs['Tissue_sub'] != 'Normal']
adata = adata[adata.obs['Tissue_sub'].cat.remove_unused_categories() != 'Blood']

adata.var=adata.var.drop(labels=['ENSEMBL', 'CHR', 'START', 'END', 'STRAND', 'BIOTYPE'],axis=1)
adata.var_names_make_unique()
adata.obs=adata.obs.drop(labels=['Barcode', 'Donor', 'Tissue', 'celltype_sub', 'Platform', 'Sample', 'Histology', 'Tissue_sub'],axis=1)

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228/adata.h5ad"
adata=sc.read_h5ad(results_file)

adata.write(results_file)

###################
#GSE138709 基因数量太少去掉
#####################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/GSE138709_RAW.tar_files')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/GSE138709_RAW.tar_files')
for x in files:
	un_gz(x)

files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/GSE138709_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/GSE138709_RAW.tar_files')
adatas = []
ss=[]
for file in files:
    adata = sc.read_csv(file)
    adata=adata.T
    strlist = file.split('_')
    ss.append(strlist[0])
    adata.obs['patients']= 'ICC' + strlist[2]
    adata.var_names_make_unique()
    adatas.append(adata)


var_names = adatas[0].var_names.intersection(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.intersection(adatas[j].var_names)

#for j in range(len(adatas)):
#    adatas[j]=re_size_adata(adatas[j],var_names)

adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=ss)
adata.obs['samples']=adata.obs.batch
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')

adata.obs=adata.obs.drop(labels=['batch'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/adata.h5ad"
adata.write(results_file)

################
#GSE125449 低于20000个基因删除
###########
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449')
files=os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449')
for x in files:
	un_gz(x)
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449/'
adatas = []
for x in ['GSE125449_Set2_','GSE125449_Set1_']:
	adata=sc.read_10x_mtx(paths,prefix=x)
	adatas.append(adata)

var_names = adatas[0].var_names.intersection(adatas[1].var_names)

#for j in range(len(adatas)):
#    adatas[j]=re_size_adata(adatas[j],var_names)

adata = adatas[0].concatenate(adatas[1],index_unique=None,batch_categories=['GSE125449_Set2_','GSE125449_Set1_'])
adata.obs['samples']=adata.obs.batch
adata.obs['patients']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
###基因太少了
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449/adata.h5ad"
adata.write(results_file)

######################################
#合并数据
######################
#2.1读取数据
#os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/LiverCancer')

adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530/adata.h5ad")
adata1.obs.samples.value_counts()

adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228/adata.h5ad")

adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449/adata.h5ad")

#sc.tl.pca(adata3)
#sc.external.pp.bbknn(adata3, batch_key='samples')
#adata3.write("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE125449/adata.h5ad")

#adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE138709/adata.h5ad")
adata1.obs_names_make_unique()
adata2.obs_names_make_unique()
adata3.obs_names_make_unique()
adata1.obs_names = adata1.obs_names + "_GSE151530"
adata2.obs_names = adata2.obs_names + "_GSE140228"
adata3.obs_names = adata3.obs_names + "_GSE125449"
#adata4.obs_names = adata4.obs_names + "_GSE138709"

duplicate_labels = adata2.obs.index.intersection(adata3.obs.index)
if len(duplicate_labels) > 0:
    print("Duplicate labels found in adata1 and adata2:", duplicate_labels)

adata = adata1.concatenate(adata2,batch_categories=["GSE151530","GSE140228"],index_unique='raise') #合并数据集

import bbknn
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.obs['datasets']=adata.obs.batch

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'LiverCancer_Integrate_cell.h5ad'
adata.write(results_file)

color_map='CMRmap'
cancer='liver'
n_comps=50
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')
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

#BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530/'
#os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/LiverCancer_Immune_cell.h5ad'
adata = sc.read_h5ad(results_file)
#adata.obs.datasets.value_counts()
#datasets
#GSE161529      71939
#E-MTAB-8107    17381
#Name: count, dtype: int64

adata = sc.read_h5ad(results_file)

#adata.obs_names = [f'{i+1}_{name}' for i, name in enumerate(adata.obs_names)]
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

import seaborn as sns
import matplotlib.pyplot as plt

#selected_classes=['B cells','T cells', 'TAMs']
#adata = adata[adata.obs.Type.isin(selected_classes)]
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)
adata.write(results_file)
```

### annotation

```python

#results_file=os.path.join(BaseDirectory, 'adata_pp_immune.h5ad')
#adata=sc.read_h5ad(results_file)
import celltypist
from celltypist import models

#models.download_models(force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"])
model_low = models.Model.load(model="Immune_All_Low.pkl")
#model_high = models.Model.load(model="Immune_All_High.pkl")
#predictions_high = celltypist.annotate(adata, model=model_high, majority_voting=True)
#predictions_high_adata = predictions_high.to_adata()
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

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_liver_celltypist_cell_label_fine')
adata.obs["louvain_1_anno"].value_counts()
for x in ['Endothelial cells','Epithelial cells']:
    adata = adata[adata.obs["louvain_1_anno"] != x]
sc.pp.pca(adata, n_comps=30)
sc.external.pp.bbknn(adata, batch_key='patients')

sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1.0)
sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

celltype_names = pivot_table.idxmax(axis=1)
adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_liver_celltypist_cell_label_fine')

adata.write(results_file)
adata=sc.read_h5ad(results_file)

```