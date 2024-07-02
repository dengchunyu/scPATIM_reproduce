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

#######################
#GSE206785
####################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE206785/')
#GSE206785_scgex.txt.gz

X = pd.read_csv('GSE206785_scgex.txt.gz',index_col=0,header=1)
n_obs=len(X.index)
n_vars = len(X.columns)
obs = pd.read_csv('GSE206785_metadata.txt.gz',sep='\t')

obs.index=X.index
var_names=X.columns
var = pd.DataFrame(index=var_names)
adata=ad.AnnData(X,obs=obs,var=var,dtype='int32')
adata.var_names_make_unique()

adata.obs['samples']= adata.obs['Sample']
adata.obs['patients']= adata.obs['Patient']
adata.obs['datasets']= 'GSE206785'

adata.obs=adata.obs.drop(labels=['Patient', 'Sample'],axis=1)

#3）输出
adata.write('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE206785/adata.h5ad')
adata.obs.Type.value_counts()
##################
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/
###############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/')
un_tar("GSE163558_RAW.tar")
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/GSE163558_RAW.tar_files/'

filesn=['GSM5004180_PT1_' ,'GSM5004181_PT2_' , 'GSM5004182_PT3_']
adatas = []
for preFix in filesn:
    adata = sc.read_10x_mtx(paths,prefix=preFix)
    adata.var_names_make_unique()
    strlist = preFix.split('_')
    adata.obs['samples']= preFix[0]
    adata.obs['patients']= strlist[1]
    adatas.append(adata)
adata = adatas[0].concatenate(adatas[1:],index_unique=None)
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='samples')
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/')
adata.write('adata.h5ad')


###################
#GSE150290
######################
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
def un_gz(file_name):
    f_name = file_name.replace(".gz", "")
    g_file = gzip.GzipFile(file_name)
    open(f_name, "wb+").write(g_file.read())
    g_file.close() #关闭gzip对象
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE150290/')
un_tar("GSE150290_RAW.tar")
files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE150290/GSE150290_RAW.tar_files/')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE150290/GSE150290_RAW.tar_files/')
files2=['GSM4546301_Pat01-B.raw_gene_bc_matrices.tar' , 'GSM4546325_Pat13-B.raw_gene_bc_matrices.tar',
'GSM4546303_Pat02-B.raw_gene_bc_matrices.tar' , 'GSM4546327_Pat14-B.raw_gene_bc_matrices.tar',
'GSM4546305_Pat03-B.raw_gene_bc_matrices.tar' , 'GSM4546329_Pat15-B.raw_gene_bc_matrices.tar',
'GSM4546307_Pat04-B.raw_gene_bc_matrices.tar' , 'GSM4546331_Pat16-B.raw_gene_bc_matrices.tar',
'GSM4546309_Pat05-B.raw_gene_bc_matrices.tar' , 'GSM4546333_Pat17-B.raw_gene_bc_matrices.tar',
'GSM4546311_Pat06-B.raw_gene_bc_matrices.tar' , 'GSM4546335_Pat18-B.raw_gene_bc_matrices.tar',
'GSM4546313_Pat07-B.raw_gene_bc_matrices.tar' , 'GSM4546337_Pat19-B.raw_gene_bc_matrices.tar',
'GSM4546315_Pat08-B.raw_gene_bc_matrices.tar' , 'GSM4546339_Pat20-B.raw_gene_bc_matrices.tar',
'GSM4546317_Pat09-B.raw_gene_bc_matrices.tar' , 'GSM4546342_Pat22-B.raw_gene_bc_matrices.tar',
'GSM4546319_Pat10-B.raw_gene_bc_matrices.tar'  ,'GSM4546344_Pat23-B.raw_gene_bc_matrices.tar',
'GSM4546321_Pat11-B.raw_gene_bc_matrices.tar' , 'GSM4546346_Pat24-B.raw_gene_bc_matrices.tar',
'GSM4546323_Pat12-B.raw_gene_bc_matrices.tar']
for i in files:
	un_gz(i)
for i in files2:
	un_tar(i)

/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE150290/GSE150290_RAW.tar_files/GSM4546301_Pat01-B.raw_gene_bc_matrices.tar_files/home/blrusti/project/01.metaplasia_scRNA/01.Cellranger/BGM-B/BGM-B/outs/raw_gene_bc_matrices/hg19/
```


### 处理细胞

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
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE206785/adata.h5ad")
for i in ["Fibroblast","Epithelial", "Endothelial","Mural","Glial"]:
	adata1 = adata1[adata1.obs["Type"] != i]
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/adata.h5ad")
adata = adata1.concatenate(adata2,batch_categories=['GSE206785','GSE163558']) #合并数据集
adata.obs_names_make_unique()
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
adata.obs['datasets']=adata.obs.batch
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
results_file = 'GastricCancer_Immune_cell.h5ad'
adata.write(results_file)

adata = sc.read_h5ad(results_file)
sc.external.pp.bbknn(adata, batch_key='patients')  

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

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)

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

def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)
adata=reshape_anno(adata)
adata.obs["louvain_1_anno"].value_counts()

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_gastric_celltypist_cell_label_fine')

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/')
adata.write(results_file)



###修改
cancer='GastricCancer'
adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/"+ cancer + "_Immune_cell.h5ad")

adata.obs.loc[adata.obs['first_celltype_annotation'].isin(['Treg']), 'first_celltype_annotation'] = 'Tcm/N.h'
marker_dotplot(cancer="GastricCancer")
sc.pl.umap(adata,color=["first_celltype_annotation"],save='_GastricCancer_first_celltype_annotation')
adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].cat.remove_unused_categories()

adata=adata[adata.obs.first_celltype_annotation != 'Mono.c']

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/"+ cancer + "_Immune_cell.h5ad")

```




