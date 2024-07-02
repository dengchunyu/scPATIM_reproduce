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
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())

from scipy.sparse import *
import numpy as np

```
```python
######################################################################################
#一、分别处理读取所有子数据集
########################
##### 1.1.E-MTAB-8107
######################

"/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer"
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/E-MTAB-8107/')
adatas = []
for i in files:
    mtx_file_counts="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/E-MTAB-8107/" + i
    adata = sc.read_csv(mtx_file_counts)
    adata = adata.T
    strlist = i.split('.')
    adata.obs['samples']= strlist[0]
    adata.obs['patients']= strlist[0]
    adata.var_names_make_unique()
    adatas.append(adata)


#var_names = adatas[0].var_names
#for j in range(len(adatas)):
#    len(adatas[j].var_names)
#    var_names = var_names.union(adatas[j].var_names)
#所有的基因数量
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

adata_concat = adata_concat[0].concatenate(adata_concat[1:],index_unique=None)
adata_concat.obs_names_make_unique()
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')

adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/E-MTAB-8107/adata.h5ad"
adata_concat.write(results_file)


###1）解压
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
#import file as f


##########################
##GSE123926 删除
##############
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926')
#1)解压缩
un_tar("GSE123926_RAW.tar")

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/GSE123926_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/GSE123926_RAW.tar_files')

for i in files:
	un_gz(i)


files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/GSE123926_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/GSE123926_RAW.tar_files')

filesn = ["GSM3516947_PDX110-","GSM3516948_PDX322-"]
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/GSE123926_RAW.tar_files'
#2）循环读取mtx文件
adatas = []
for preFix in filesn:
	adata = sc.read_10x_mtx(paths,prefix=preFix)
	adata.var_names_make_unique()
	adatas.append(adata)

var_names = adatas[0].var_names.union(adatas[1].var_names)

for j in range(39):
    len(adatas[j].var_names)
    var_names = var_names.union(adatas[j].var_names)


#3）整合以上样本
adata = adatas[0].concatenate(adatas[1:],index_unique=None)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/adata.h5ad"
adata.write(results_file)

###################
##GSE161529
##############
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529')
#1)解压缩
import tarfile
tar = tarfile.open("GSE161529_RAW.tar")
tar.extractall()
tar.close()
#un_tar("GSE161529_RAW.tar")


files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/GSE161529_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/GSE161529_RAW.tar_files')
####R
files_pre<-read.table("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/file_name.txt")
a<-unique(files_pre$V1)
a<-paste0(a,"-")
a<-data.frame(a)
write.table(a,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/file_name.txt",quote=F,row.names=F,col.names=F)
######
files=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/file_name.txt")
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/GSE161529_RAW.tar_files'
#2）循环读取mtx文件
import shutil
from shutil import copyfile
filesn = ['GSM4909317_ER-MH0173-T-', 'GSM4909312_ER-MH0056-LN-', 'GSM4909316_ER-MH0167-LN-', 'GSM4909292_HER2-MH0069-', 'GSM4909290_HER2-PM0337-', 'GSM4909285_TN-B1-MH4031-', 'GSM4909311_ER-MH0056-T-',  'GSM4909288_TN-B1-MH0177-', 'GSM4909306_ER-MH0029-9C-',  'GSM4909308_ER-MH0040-LN-', 'GSM4909293_HER2-MH0161-', 'GSM4909300_ER-MH0032-', 'GSM4909318_ER-MH0173-LN-', 'GSM4909310_ER-MH0043-LN-', 'GSM4909314_ER-MH0064-LN-', 'GSM4909297_ER-MH0125-', 'GSM4909299_ER-MH0114-T3-', 'GSM4909291_HER2-MH0031-', 'GSM4909281_TN-MH0126-', 'GSM4909319_mER-PM0178-', 'GSM4909298_ER-PM0360-', 'GSM4909282_TN-MH0135-', 'GSM4909294_HER2-MH0176-', 'GSM4909320_mER-MH0068-T-', 'GSM4909283_TN-SH0106-', 'GSM4909296_ER-MH0001-', 'GSM4909302_ER-MH0025-', 'GSM4909309_ER-MH0043-T-', 'GSM4909303_ER-MH0151-', 'GSM4909304_ER-MH0163-', 'GSM4909305_ER-MH0029-7C-', 'GSM4909301_ER-MH0042-', 'GSM4909287_TN-B1-Tum0554-', 'GSM4909286_TN-B1-MH0131-', 'GSM4909307_ER-MH0040-', 'GSM4909295_ER-AH0319-', 'GSM4909313_ER-MH0064-T-', 'GSM4909284_TN-MH0114-T2-', 'GSM4909321_mER-MH0068-LN-', 'GSM4909289_HER2-AH0308-']
adatas = []
for preFix in filesn:
    copyfile("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/GSE161529_features.tsv.gz", "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/GSE161529_RAW.tar_files/"+preFix+"features.tsv.gz")

for preFix in filesn:
    adata = sc.read_10x_mtx(paths,prefix=preFix)
    adata.var_names_make_unique()
    strlist = preFix.split('_')
    adata.obs['samples']= preFix
    adata.obs['patients']= strlist[1].split('-')[2]
    adatas.append(adata)

#3）整合以上样本
adata = adatas[0].concatenate(adatas[1:],index_unique=None)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='samples')

results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/adata.h5ad"
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
adata.var=adata.var.drop(labels=['feature_types'],axis=1)
tumor_data.obs.sample_id = tumor_data.obs.samples.str.split('_').str[0]
tumor_data.obs.tumor_type = tumor_data.obs.samples.str.split('_').str[1].str.split('-').str[0]
tumor_data.obs.patient = tumor_data.obs.samples.str.split('_').str[1].str.split('-').str[1]
tumor_data.obs.tumor_type.value_counts()
adata.write(results_file)

###################
##GSE176078 重点，有注释的数据
##############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078')
#1)解压缩
un_gz("GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")
un_tar("GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar")

paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar_files/Wu_etal_2021_BRCA_scRNASeq/'
os.chdir(paths)

adata = sc.read_mtx('count_matrix_sparse.mtx')
adata=adata.T
genes = pd.read_csv('count_matrix_genes.tsv',names=['genes'])
adata.var_names=genes['genes']
barcodes = pd.read_csv('count_matrix_barcodes.tsv',names=['barcodes'])
adata.obs_names=barcodes['barcodes']
metadata = pd.read_csv('metadata.csv', index_col=0)
adata.obs = metadata
adata = adata[adata.obs_names, adata.var_names]
adata.obs['samples']= adata.obs['orig.ident']
adata.obs['patients']= adata.obs['orig.ident']
adata.obs=adata.obs.drop(labels=['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mito', 'subtype', 'celltype_subset', 'celltype_minor',],axis=1)

#3）输出
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/adata.h5ad"
adata.write(results_file)

######################################
#合并数据
######################
#2.1读取数据
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')

#adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/adata.h5ad")
adata1.obs['samples']= adata1.obs['orig.ident']
adata1.obs['patients']= adata1.obs['orig.ident']
adata1.obs=adata1.obs.drop(labels=['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mito', 'subtype', 'celltype_subset', 'celltype_minor',],axis=1)
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/adata.h5ad")
#adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE123926/adata.h5ad")
adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/E-MTAB-8107/adata.h5ad")
adata4.obs_names_make_unique()
adata2.obs_names_make_unique()

adata_concat = adata1.concatenate(adata2,adata4, batch_categories=['GSE176078','GSE161529','E-MTAB-8107']) #合并数据集

sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
adata_concat.obs['datasets']=adata_concat.obs.batch
sc.tl.pca(adata_concat, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='breast'
n_comps=50
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'breastcancer')

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'breastcancer_integrate_cell.h5ad'
adata_concat.write(results_file)

import scipy
from scipy.sparse import *
import numpy as np

#isspmatrix_csc(x)
########################
#2.3 给予adata1作为注释数据，将其他数据的细胞类型初步注释出来
var_names=adata1.var_names.intersection(adata2.var_names)
var_names=var_names.intersection(adata4.var_names)
adata1=adata1[:,var_names]
adata2=adata2[:,var_names]
adata4=adata4[:,var_names]

sc.pp.pca(adata1)
sc.pp.neighbors(adata1)
sc.tl.umap(adata1)

sc.tl.ingest(adata2, adata1, obs='celltype_major')
#sc.tl.ingest(adata3, adata1, obs='celltype_major')
sc.tl.ingest(adata4, adata1, obs='celltype_major')

adata_concat = adata1.concatenate(adata2,adata4, batch_categories=['GSE176078','GSE161529','E-MTAB-8107']) #合并数据集
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
adata_concat.obs['datasets']=adata_concat.obs.batch

adata_concat.obs=adata_concat.obs.drop(labels=['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mito', 'subtype', 'celltype_subset', 'celltype_minor',],axis=1)

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'breastcancer_integrate_cell.h5ad'
adata_concat.write(results_file)

#############################
#2.5 删除无用细胞
adata_concat = adata_concat[adata_concat.obs.celltype_major != 'Normal Epithelial']
adata_concat = adata_concat[adata_concat.obs.celltype_major != 'Cancer Epithelial']
adata_concat = adata_concat[adata_concat.obs.celltype_major != 'CAFs']
adata_concat = adata_concat[adata_concat.obs.celltype_major != 'Endothelial']
adata_concat = adata_concat[adata_concat.obs.celltype_major != 'PVL']
adata_concat.obs=adata_concat.obs.drop(labels=['celltype_major','batch'],axis=1)

#os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/Pancancer/Immune_celldata')
#############################
adata_concat.write('BreastCancer_Immune_cell.h5ad')

adata_concat = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/BreastCancer_Immune_cell.h5ad')
#adata_concat.obs_names_make_unique()
#sc.pp.normalize_total(adata_concat, target_sum=1e4) ##标准化
#sc.pp.log1p(adata_concat)
#adata_concat.write('breastcancer_Immune_cell.h5ad')
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
import scarches as sca
import urllib.request

BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/'

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/BreastCancer_Immune_cell.h5ad")

os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='BreastCancer_Immune_cell.h5ad'

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

sc.pp.filter_genes(adata, min_cells=50)
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=4000)
sc.pp.filter_cells(adata, max_counts=50000)

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
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1.0)
class_counts = adata.obs.louvain_res1.value_counts()
selected_classes = class_counts[class_counts >= 100].index.tolist()
adata = adata[adata.obs.louvain_res1.isin(selected_classes)]

```

### celltypist

```python
#BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/'
#os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

model_low = models.Model.load(model="Immune_All_Low.pkl")
duplicated_obs_names = adata.obs_names.duplicated(keep='first')
adata = adata[~duplicated_obs_names]

adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

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


sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_breast_celltypist_cell_label_fine')

for x in ['Late erythroid','Epithelial cells']:
    adata = adata[adata.obs["louvain_1_anno"] != x]
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1.0)
sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]

pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)

celltype_names = pivot_table.idxmax(axis=1)
adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_breast_celltypist_cell_label_fine')

adata.write(results_file)
##########################################


```


old code 

```python

marker_genes = {"Tn":["LEF1", "TCF7", "KLF2", "TXK", "BACH2","LTB", "FLT3LG", "TNFSF8", "CMTM8", "IL23A", "TIMP1", "WNT7A","CCR7", "IL7R", "IL6R", "IFNGR2","SELL", "MAL", "EEF1A1", "ACTN1", "TRABD2A", "TPT1", "EEF1B2", "NELL2", "NOSIP", "PABPC1"],
                "GZMK+ Tem":["EOMES", "SUB1", "GZMK", "GZMA", "CCL5", "GZMH", "IL32", "CCL4", "CD74", "CCR5", "CXCR3","HLA-DRB1", "HLA-DPA1", "COTL1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB5", "ITM2C", "CST7", "APOBEC3G"],
                "Temra":["ASCL2", "KLF2", "KLF3", "ZEB2","TBX21", "GZMH", "GZMB", "GZMM", "GZMA", "CX3CR1", "CXCR2", "CMKLR1", "CXCR1", "FGFBP2", "FCGR3A", "S1PR5", "PRSS23", "GNLY", "NKG7","KLRD1", "FGR", "PLEK", "C1orf21"],
                "ZNF683+CXCR6+ Trm":["ZNF683", "HOPX", "ID2","ZFP36L2", "RBPJ", "CKLF", "IL32", "GZMB", "XCL1", "CCL5", "GZMA","XCL2", "CXCR6", "CXCR3", "CAPG", "TMSB4X", "S100A4", "LGALS3", "ACTB", "SH3BGRL3", "CD52","LGALS1", "ITGA1", "LDLRAD4"],
                "Mast": "TPSAB1",
                "GZMK+ Tex":["TOX", "TSC22D1", "PRDM1", "TRPS1", "EOMES", "GZMK", "CXCL13", "GZMA", "CCL3", "CCL5", "CCL3L3", "TNFSF4", "CCL4", "IFNG", "NAMPT", "CD74", "CXCR6", "CCR5", "HAVCR2", "CD27", "VCAM1", "LYST", "PDCD1", "DUSP4", "CTLA4", "TNFRSF9", "HLA-DQA1", "HLA-DRB1"],
                "terminal Tex":["RBPJ", "ETV1", "TOX", "ZBED2", "TOX2", "CXCL13", "TNFSF4", "FAM3C", "GZMB", "CSF1", "CCL3", "CD70", "IFNG", "NAMPT", "FASLG", "IL2RA", "CXCR6", "CD74", "IL2RB", "IL2RG", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4", "KRT86", "TNFRSF18", "GEM", "TIGIT", "DUSP4"]
                }
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

selected_classes=['1','2','4','5','6','7','10','12']

adataT = adata[adata.obs.louvain_res1.isin(selected_classes)]

sc.pl.dotplot(adataT,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('T2_marker_dotplot.pdf', format="pdf")
plt.clf()
marker_genes = {'TN':['CCR7','PABPC1','KLF2','IL7R','PASK','FTH1','SELL','S1PR1','EEF1B2','RPS8','RPS12','RPL3','RPL32','RPS3A','CD55','RPL34','RPS13','EEF1G','RPL13','EEF1A1','RPS18','GPR183','RPS6','RPL11','RPL5','GIMAP7','SLC2A3','RPL22','MT1X','MT2A'],
'TH17':['KLRB1' ,'TRAV1-2' ,'RORA' ,'ERN1' ,'CTSH' ,'IL4I1' ,'PERP' ,'CD40LG' ,'AQP3' ,'MGAT4A' ,'IL7R' ,'FKBP11' ,'CCL20' ,'TMIGD2' ,'SYTL2' ,'S100A6' ,'ADAM19' ,'FURIN' ,'OSTF1' ,'S100A4' ,'TNFRSF25' ,'JAML' ,'ALOX5AP' ,'LGALS3' ,'CXCR6' ,'PIM1' ,'CAPG' ,'LTB' ,'CD69' ,'PLIN2'],
'TFH1':['CXCL13' ,'KLRB1' ,'NR3C1' ,'CTSL' ,'CTSH' ,'PTPN13' ,'TNFSF8' ,'CD40LG' ,'UNQ6494' ,'CD200' ,'KIAA0319L' ,'GNA15' ,'RBPJ' ,'SPOCK2' ,'ADAM19' ,'MAF' ,'FKBP5' ,'CD6' ,'GK' ,'SNX9' ,'TMEM173' ,'CCDC50' ,'TNFRSF25' ,'FURIN' ,'RORA' ,'BHLHE40' ,'NMB' ,'PLIN2' ,'CD83' ,'CCL20'],
'TFH2':['NMB' ,'TOX2' ,'PASK' ,'NR3C1' ,'GNG4' ,'IL6ST' ,'SESN3' ,'PTPN13' ,'FKBP5' ,'TCF7' ,'PGM2L1' ,'FAAH2' ,'GK' ,'CD200' ,'MAGEH1' ,'IGFBP4' ,'CXCR5' ,'TNFSF8' ,'SMCO4' ,'ICA1' ,'TSHZ2' ,'CHI3L2' ,'TOX' ,'CNIH1' ,'RNF19A' ,'SH2D1A' ,'PPP1CC' ,'WDR74' ,'CXCL13' ,'JUNB'],
'TMEM-CD4':['LMNA' ,'ANXA1' ,'IL7R' ,'CCR7' ,'FOS' ,'ZFP36' ,'GPR183' ,'MYADM' ,'CDKN1A' ,'FTH1' ,'VIM' ,'S100A10' ,'SLC2A3' ,'ZFP36L2' ,'CXCR4' ,'PABPC1' ,'TPT1' ,'S1PR1' ,'RPL3' ,'KLF2' ,'FOSB' ,'AHNAK' ,'RGCC' ,'CD55' ,'PTGER4' ,'TUBB2A' ,'KLF6' ,'NFKBIA' ,'SOD2' ,'MT1X'],
'Treg':['TNFRSF4' ,'FOXP3' ,'IL2RA' ,'LAIR2' ,'BATF' ,'TNFRSF18' ,'LTB' ,'IL1R2' ,'SAT1' ,'CCR8' ,'IL32' ,'CTLA4' ,'MIR4435-2HG' ,'GK' ,'ICOS' ,'AC133644.2' ,'MAGEH1' ,'S100A4' ,'TBC1D4' ,'CARD16' ,'HTATIP2' ,'DNPH1' ,'PIM2' ,'CORO1B' ,'CTSC' ,'AC145110.1' ,'GLRX' ,'SYNGR2' ,'CCR6' ,'TNFRSF1B'],
'TEFF':['GZMK' ,'CCL4L2' ,'TUBA4A' ,'CST7' ,'GZMH' ,'ITM2C' ,'PIK3R1' ,'DUSP2' ,'CCL4' ,'NKG7' ,'ANXA1' ,'LITAF' ,'KLRG1' ,'RP11-291B21.2' ,'CRTAM' ,'CCL5' ,'CXCR4' ,'GZMM' ,'XCL2' ,'AOAH' ,'MYADM' ,'CD8A' ,'CD8B' ,'TC2N' ,'NFATC2' ,'PLEK' ,'HLA-DPB1' ,'CMC1' ,'MT1F' ,'MT1X'],
'TMEM-CD8':['ANXA1' ,'LMNA' ,'MYADM' ,'ZFP36L2' ,'IL7R' ,'S100A10' ,'ZFP36' ,'VIM' ,'CD8A' ,'CD8B' ,'ZNF683' ,'CCL5' ,'CXCR4' ,'PTGER4' ,'PARP8' ,'CD55' ,'S100A6' ,'TNFAIP3' ,'AHNAK' ,'XCL1' ,'TUBB2A' ,'GPR65' ,'BTG2' ,'CD69' ,'GPR183' ,'RGCC' ,'FOSB' ,'TOB1' ,'FOS' ,'LDLRAD4'],
'TEX':['CXCL13' ,'CCL3' ,'CD8A' ,'GZMB' ,'VCAM1' ,'NKG7' ,'CCL5' ,'HLA-DRA' ,'KRT86' ,'LAG3' ,'CCL4' ,'RGS2' ,'RHOB' ,'CD8B' ,'HLA-DQA1' ,'IFNG' ,'HAVCR2' ,'RP11-291B21.2' ,'GZMH' ,'HLA-DRB1' ,'STMN1' ,'GNLY' ,'PLPP1' ,'HLA-DRB5' ,'ALOX5AP' ,'GZMA' ,'FABP5' ,'LRRN3' ,'CCL4L2' ,'CTSW'],
'NK/NKT':['TYROBP' ,'FCER1G' ,'XCL1' ,'AREG' ,'GNLY' ,'XCL2' ,'TRDC' ,'KLRC1' ,'GZMA' ,'SH2D1B' ,'CTSW' ,'KLRD1' ,'IFITM3' ,'KRT81' ,'HOPX' ,'LAT2' ,'MATK' ,'ADGRG3' ,'NR4A1' ,'TMIGD2' ,'KIR2DL4' ,'KRT86' ,'GSTP1' ,'CD83' ,'TXK' ,'KLRB1' ,'B3GNT7' ,'GADD45B' ,'CD247' ,'IL2RB']
}
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found


sc.pl.dotplot(adataT,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('T_marker_dotplot.pdf', format="pdf")
plt.clf()

marker_genes = { "Temra":["ASCL2", "KLF2", "KLF3", "ZEB2", "TBX21","GZMH", "GZMB", "GZMM", "GZMA","CX3CR1", "CXCR2", "CMKLR1", "CXCR1","FGFBP2", "FCGR3A", "S1PR5", "PRSS23", "GNLY", "NKG7", "KLRD1","FGR", "PLEK", "C1orf21"],
                "CD4+ T naive": ["CD4", "IL7R", "TRBC2", "CCR7"],
                "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
                "T activation": ["CD69", "CD38"],
                "T naive": ["LEF1", "CCR7", "TCF7"],
                "G/M prog": ["MPO", "BCL2", "KCNQ5", "CSF3R"],
                "HSC": ["NRIP1", "MECOM", "PROM1", "NKAIN2", "CD34"]}
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

sc.pl.dotplot(adataT,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('T3_marker_dotplot.pdf', format="pdf")
plt.clf()

method="wilcoxon"
sc.tl.rank_genes_groups(adataT, groupby="louvain_res1", method="wilcoxon", key_added="diff_annotation")

sc.tl.filter_rank_genes_groups(adataT,min_in_group_fraction=0.2,max_out_group_fraction=0.2, key="diff_annotation",key_added="diff_annotation_filtered")

sc.tl.dendrogram(adataT,groupby="louvain_res1")
sc.pl.rank_genes_groups_dotplot(adataT,groupby="louvain_res1",standard_scale="var",n_genes=5,key="diff_annotation_filtered",save='diff_T_louvain_res1_filtered_dotplot')
###memorycell
m_marker={"TM":["ZFP36L2", "KLF2", "KLF3", "TCF7", "HOPX", "LMO4"]}

sc.pl.dotplot(adataT,groupby="louvain_res1",var_names=m_marker,standard_scale="var")
plt.savefig('T4_marker_dotplot.pdf', format="pdf")
plt.clf()



marker_genes = {"cDC1": ["CLEC9A", "CADM1","TNFRSF1B","HMGN2","SMAP2","RHOC","PSMD4","ETV3","CD1E","RCSD1","MPZL1","GLUL","RGS1","G0S2","OPN3","RHOB","PPM1B","REL","PCBP1","MOB1A","CAPG","MERTK","FOXD4L1","CXCR4","ZEB2","TFPI","RAPH1","FILIP1L","NFKBIZ","ALCAM","CD86","TPRA1","FNDC3B","CCDC50","RHOH","ANKRD55","CDC42SE2","IL12B","MGAT1","LINC00847","DUSP22","FLOT1","IER3","PHIP","TNFAIP3","PLEKHG1"],
                "cDC2": ["CST3","COTL1","LYZ","DMXL2","CLEC10A","FCER1A","TNFRSF14","TXNIP","SLAMF1","SLAMF7","F11R","TSTD1","MGST3","TNFSF4","TOR3A","RAB29","HLX","GPR137B","VAMP5","IGKV1-5","IGKV1-17","IGKV3-20","BCL2L11","SLC20A1","GYPC","TNFAIP6","IFIH1","WIPF1","CALCRL","CNOT10","H1FX","SLC9A9","LAP3","CXCL11","SEPT11","RP11-290F5.1","LMNB1","HINT1","TIFAB","ADAM19","CREBRF","ADTRP","HLA-DOB","TAP1","TAPBP","BZW2","IKZF1","SRI","TSPAN33","PIM2","PPP3CC","DPYSL2","RAB11FIP1","RCL1","DAPK1","PTBP3","LSP1","PACS1","GSTP1","SERPINH1","FAM21C","PPA1","COMTD1","PTMS","C12orf57","ETV6","SELPLG","RHOF","ATP6V0A2","FLT3"], 
                "pDC": ["IL3RA", "COBLL1", "TCF4", "PTGDS","SOX4","GZMB","IRF7", "LILRA4", "DERL3","TPM2", "JCHAIN","MZB1", "PTCRA"],
                "M1": ["IL23","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5","IRF5","IRF1","CD40","IDO1","KYNU","CCR7"],
                "M2":["IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1","VEGFA","VEGFB","VEGFC","VEGFD","EGF","CTSA","CTSB","CTSC","CTSD","TGFB1","TGFB2","TGFB3","MMP14","MMP19","MMP9","CLEC7A","WNT7B","FASL","TNFSF12","TNFSF8","CD276","VTCN1","MSR1","FN1","IRF4"],
                "Angiogenesis_mac":["CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA"]
}

marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

selected_classes=['0','8','9','11','13','14','15','16','18']

adataM = adata[adata.obs.louvain_res1.isin(selected_classes)]

sc.pl.dotplot(adataM,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('M1_marker_dotplot.pdf', format="pdf")
plt.clf()

method="wilcoxon"
sc.tl.rank_genes_groups(adataM, groupby="louvain_res1", method="wilcoxon", key_added="diff_annotation")

sc.tl.filter_rank_genes_groups(adataM,min_in_group_fraction=0.2,max_out_group_fraction=0.2, key="diff_annotation",key_added="diff_annotation_filtered")

sc.tl.dendrogram(adataM,groupby="louvain_res1")
sc.pl.rank_genes_groups_dotplot(adataM,groupby="louvain_res1",standard_scale="var",n_genes=5,key="diff_annotation_filtered",save='diff_M_louvain_res1_filtered_dotplot')

marker_genes = {'MONOCYTE': ['C1QC','FCN1','S100A6','C1QB','C1QA','APOE','TREM2','S100A4','APOC1','S100A9','LYZ','A2M','S100A8','VCAN','SH3BGRL3','CD52','GPR34','TIMP1','C3','CD74','H3F3A','EMP3','PLTP','ITM2B','AKR1B1','CD59','OLFML3','RPS9','LST1','STXBP2','RPL28','S100A12','FCGBP','S100A10','SLCO2B1','RP11-1143G9.4','CFP','CSTA','BIN1','FLNA','RPL39','MSR1','VSIG4','CD300E','AXL','LGALS3BP','BHLHE41','LINC01272','COTL1']}
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

sc.pl.dotplot(adataM,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('M7_marker_dotplot.pdf', format="pdf")
plt.clf()

###
#T注释结果
cl_annotation = {"0": "Monocyte",
    "1": "Tnaive/mem",
    "2": "Tnaive/mem",
    "3":"Bcells",
    "4": "NKT",
     "5": "Teff",
     "6":"Treg",
     "7":"Tex",
     "8":"cDC",
     "9":"IFN.TAM",
     "10":"Temra",
     "11":"FCGBP+TAM",
     "12":"Treg",
     "13":"CXCL9+TAM",
     "14":"CCL8+TAM",
     "15":"Hypoxic.TAM",
     "16":"pDC",
     "17":"Plasma",
    "18":"prol.TAM"}
adata.obs["manual_celltype_annotation"] = adata.obs.louvain_res1.map(cl_annotation)

adata.write(results_file)

method="wilcoxon"
sc.tl.rank_genes_groups(adata, groupby="manual_celltype_annotation", method="wilcoxon", key_added="dea_1_annotation")
sc.tl.filter_rank_genes_groups(adata,min_in_group_fraction=0.1,max_out_group_fraction =0.5, key="dea_1_annotation",key_added="dea_1_annotation_filtered")


n_genes=500
cluster_method="wilcoxon"
result = adata.uns['dea_1_annotation']
groups = result['names'].dtype.names
pval_table = pd.DataFrame({group + '_' + key[:2]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join(BaseDirectory, method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_1_annotation.xlsx'), engine='openpyxl')

result = adata.uns['dea_1_annotation_filtered']
groups = result['names'].dtype.names
pval_table = pd.DataFrame({group + '_' + key[:2]: result[key][group] for group in groups for key in ['names', 'pvals_adj']}).head(n_genes)
pval_table.to_excel(os.path.join(BaseDirectory, method + '_pval_table_' + cluster_method + '_clusters_' + str(n_genes) + 'genes_1_annotation_filtered.xlsx'), engine='openpyxl')

sc.tl.dendrogram(adata,groupby="manual_celltype_annotation")
sc.pl.rank_genes_groups_dotplot(adata,groupby="manual_celltype_annotation",standard_scale="var",n_genes=5,key="dea_1_annotation_filtered",save='diff_manual_celltype_filtered_dotplot')

adata.write(results_file)

sc.pl.umap(adata,color=["manual_celltype_annotation"],save='_louvain_annotation')

```