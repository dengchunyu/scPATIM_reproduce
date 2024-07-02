```python
import tarfile
import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import numpy as np
import warnings
import gc
from anndata import AnnData
import matplotlib.pyplot as pt
from scipy.sparse import *
import anndata as ad
import scipy
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())
# load adata

######################################################################################
#一、分别处理读取所有子数据集
########################
##### 1.1.GSE155698
######################

os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698')
import gzip

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


un_tar("GSE155698_RAW.tar")

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files')

for i in files:
	un_gz(i)

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files')

for i in files:
	un_tar(i)

# 2）移动文件夹
import shutil

def remove_file(old_path, new_path):
    print(old_path)
    print(new_path)
    filelist = os.listdir(old_path) #列出该目录下的所有文件,listdir返回的文件列表是不包含路径的。
    print(filelist)
    for file in filelist:
        src = os.path.join(old_path, file)
        dst = os.path.join(new_path, file)
        print('src:', src)
        print('dst:', dst)
        shutil.move(src, dst)

# 数据前移       
path = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files'
for rootP in (os.listdir(path)):
	rootPath = os.path.join(path,rootP)
	filesList = os.listdir(rootPath)
	for file in filesList:
		filePath=os.path.join(rootPath,file)
		if os.path.isdir(filePath):
			remove_file(filePath,rootPath)

###这里因为数据包含类型比较杂乱，有h5文件类型以及mtx文件类型
#3）读取h5 files
pdac_14 = sc.read_10x_h5(filename="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files/GSM4710703_PDAC_TISSUE_14.tar_files/filtered_feature_bc_matrix.h5", genome=None, gex_only=True, backup_url=None)
results_file = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/pdac_14.h5ad' ##置结果文件保存路径
pdac_14.write(results_file)

pdac_16 = sc.read_10x_h5(filename="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files/GSM4710705_PDAC_TISSUE_16.tar_files/filtered_gene_bc_matrices_h5.h5", genome=None, gex_only=True, backup_url=None)
results_file = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/pdac_16.h5ad' ##置结果文件保存路径
pdac_16.write(results_file)

#4.）读取mtx files
path = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files'
filelist = os.listdir(path)
adatas = []
ss=[]
pp=[]
for i in filelist:
    adata = sc.read_10x_mtx(path="/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/GSE155698_RAW.tar_files/" + i +"/filtered_feature_bc_matrix/")
    strlist = i.split('.')
    strlist2 = strlist[0].split('_')
    ss.append(strlist2[0])
    pp.append(strlist2[3])
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.obs['patients']=strlist2[3]
    adatas.append(adata)

adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/pdac_14.h5ad")
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.obs['patients']='14'
adatas.append(adata)
ss.append('GSM4710703')

adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/pdac_16.h5ad")
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.obs['patients']='16'
adatas.append(adata)
ss.append('GSM4710705')



#5）对以上的多样本进行整合
adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=ss)
adata.obs['samples']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)

adata.var=adata.var.drop(labels=['feature_types-GSM4710689', 'feature_types-GSM4710690', 'feature_types-GSM4710691', 'feature_types-GSM4710692', 'feature_types-GSM4710693', 'feature_types-GSM4710694', 'feature_types-GSM4710695', 'feature_types-GSM4710696', 'feature_types-GSM4710697', 'feature_types-GSM4710698', 'feature_types-GSM4710699', 'feature_types-GSM4710700', 'feature_types-GSM4710701', 'feature_types-GSM4710702', 'feature_types-GSM4710703', 'genome-GSM4710703', 'feature_types-GSM4710704'],axis=1)

#6）输出数据
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE155698/adata.h5ad"
adata.write(results_file)

####################################
##### 1.2.GSE154778
####################################
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778')
#1)解压缩
un_tar("GSE154778_RAW.tar")

files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/GSE154778_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/GSE154778_RAW.tar_files')

for i in files:
	un_gz(i)


files = os.listdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/GSE154778_RAW.tar_files')
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/GSE154778_RAW.tar_files')

filesn = ["GSM4679532_K16733_","GSM4679533_Y00006_","GSM4679534_T2_","GSM4679535_T3_","GSM4679536_T4_","GSM4679537_T5_","GSM4679538_T6_","GSM4679539_T8_","GSM4679540_T9_","GSM4679541_T10_"]
paths='/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/GSE154778_RAW.tar_files'
#2）循环读取mtx文件
adatas = []
ss=[]
for preFix in filesn:
    adata = sc.read_10x_mtx(paths,prefix=preFix)
    strlist = preFix.split('_')
    ss.append(strlist[0])
    adata.obs['patients']=strlist[1]
    adata.var_names_make_unique()
    adatas.append(adata)

var_names = adatas[0].var_names.union(adatas[1].var_names)
for j in range(len(adatas)):
    len(adatas[j].var_names)
    var_names = var_names.union(adatas[j].var_names)

for j in range(len(adatas)):
    adatas[j]=re_size_adata(adatas[j],var_names)

#3）整合以上样本
adata = adatas[0].concatenate(adatas[1:],index_unique=None,batch_categories=ss)
adata.obs['samples']=adata.obs.batch
adata.obs=adata.obs.drop(labels=['batch'],axis=1)
results_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/GSE154778/adata.h5ad"
adata.write(results_file)

#####################################
## 1.3.SCP1644
####################################
#SCP1644只有一种文件类型
mtx_file_counts = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/SCP1644/Biopsy_RawDGE_23042cells.csv"
sample_metadata_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/SCP1644/complete_MetaData_70170cells_scp.csv"
out_file = "/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/SCP1644/adata.h5ad"
#1）读取csv文件矩阵
print('Reading metadata')
metadata = pd.read_csv(sample_metadata_file, index_col=0)
# Read mtx file as anndata
print('Reading count matrix')
dataset_counts = sc.read_csv(mtx_file_counts)
print('Transposing matrices')
dataset_counts = dataset_counts.T
# Append metadata
print('Appending metadata')
#2）整合meta数据与矩阵
dataset_counts.obs = metadata.loc[dataset_counts.obs_names]
dataset_counts = dataset_counts[dataset_counts.obs_names, dataset_counts.var_names]
dataset_counts.obs['samples']=dataset_counts.obs.biosample_id
dataset_counts.obs['patients']=dataset_counts.obs.donor_ID

dataset_counts.obs.Coarse_Cell_Annotations

dataset_counts.obs=dataset_counts.obs.drop(labels=['biosample_id', 'donor_ID', 'Genes', 'UMI', 'dataset', 'sample.type', 'exvivo.treatment', 'species', 'species__ontology_label', 'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label', 'sex'],axis=1)

#3）输出
dataset_counts.write(out_file)


########################################################################################
#二、合并所有的子数据集
#########################################################################################
def re_size_adata(adata,var_names):
    var_name=var_names.difference(adata.var_names)
    len(var_name)
    n1=len(var_name)
    n2=len(adata.obs_names)
    z=np.zeros((n2,n1))
    z=csc_matrix(z,dtype='float64')
    if isspmatrix_csc(adata.X):
        X = scipy.sparse.hstack((adata.X,z))
    else:
        X = csc_matrix(adata.X,dtype='float64')
        X = hstack((X,z))
    var_names= adata.var_names.union(var_name)
    var = pd.DataFrame(index=var_names)
    adata=ad.AnnData(X,obs=adata.obs,var=var,dtype='float64')
    return adata

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
######################
#2.1读取数据
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/SCP1644/adata.h5ad")
sc.tl.pca(adata1)
sc.external.pp.bbknn(adata1, batch_key='patients')
#adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE154763/adata.h5ad")
#adata2=adata2[adata2.obs['tissue']=='T']
#adata2.obs['patients']=adata2.obs['patient']
#adata2.obs['samples']=adata2.obs['patient']
#adata2.obs=adata2.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch','library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster','source', 'tech', 'UMAP1', 'UMAP2'],axis=1)

adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE155698/adata.h5ad")
ps=[]
for x in adata3.obs.patients:
    x="PancreasPatient_"+str(x)
    ps.append(x)

adata3.obs.patients= ps
adata3.obs.samples= ps
sc.tl.pca(adata3)
sc.external.pp.bbknn(adata3, batch_key='patients')

adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE154778/adata.h5ad")
sc.tl.pca(adata4)
sc.external.pp.bbknn(adata4, batch_key='patients')

adata5 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE156728/adata.h5ad")

adata5.obs['patients']=adata5.obs['patient']
adata5.obs['samples']=adata5.obs['patient']
adata5.obs=adata5.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)
sc.tl.pca(adata5)
sc.external.pp.bbknn(adata5, batch_key='patients')

#####################
#2.2 得到所有子数据集基因交集
adata3.obs_names_make_unique()
adata4.obs_names_make_unique()
adata5.obs_names_make_unique()

adata_concat = adata1.concatenate(adata3,adata4,adata5, batch_categories=['SCP1644', 'GSE155698','GSE154778','GSE156728']) #合并数据集
import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
adata_concat.obs['datasets']=adata_concat.obs.batch


os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Pancreas_integrate_cell.h5ad'
adata_concat.write(results_file)

adata_concat.obs.Coarse_Cell_Annotations.unique()[10:15]


for i in ['Tumor','Endothelial','Mesenchymal', 'Hepatocyte']:
    adata_concat = adata_concat[adata_concat.obs.Coarse_Cell_Annotations != i]

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Pancreas_Immune_cell.h5ad'

adata_concat.obs=adata_concat.obs.drop(labels=['Coarse_Cell_Annotations'],axis=1)
adata_concat.obs['datasets']=adata_concat.obs['batch']
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)
adata_concat.write(results_file)

sc.tl.pca(adata_concat, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='Pancreatic'
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_'+ 'comps_PCA')

#var_names = adata1.var_names.intersection(adata2.var_names) 

var_names = adata1.var_names.union(adata3.var_names)
#var_names = var_names.union(adata3.var_names)
var_names = var_names.union(adata4.var_names)
var_names = var_names.union(adata5.var_names)

adata1=re_size_adata(adata1,var_names)
#adata2=re_size_adata(adata2,var_names)
adata3=re_size_adata(adata3,var_names)
adata4=re_size_adata(adata4,var_names)
adata5=re_size_adata(adata5,var_names)


#adata1 = adata_concat[adata_concat.obs.batch == 'SCP1644']
########################
#2.3 给予adata1作为注释数据，将其他数据的细胞类型注释出来


sc.pp.pca(adata1)
sc.pp.neighbors(adata1)
sc.tl.umap(adata1)
#sc.pl.umap(adata1, color='Coarse_Cell_Annotations',legend_loc='on data',legend_fontsize='xx-small',save='pancreas_SCP1644_umap.pdf')

#sc.tl.ingest(adata2, adata1, obs='Coarse_Cell_Annotations')
sc.tl.ingest(adata3, adata1, obs='Coarse_Cell_Annotations')
sc.tl.ingest(adata4, adata1, obs='Coarse_Cell_Annotations')
sc.tl.ingest(adata5, adata1, obs='Coarse_Cell_Annotations')

#############################
#2.4 整合这些子数据集

adata_concat = adata1.concatenate(adata3,adata4,adata5, batch_categories=['SCP1644', 'GSE155698','GSE154778','GSE156728']) #合并数据集
import bbknn
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
#############################
#2.6 输出结果
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Pancreas_integrate_cell.h5ad'
adata_concat.write(results_file)

adata_concat.obs.Coarse_Cell_Annotations.unique()[10:15]

for i in ['Tumor','Endothelial','Mesenchymal', 'Hepatocyte']:
    adata_concat = adata_concat[adata_concat.obs.Coarse_Cell_Annotations != i]

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell')
results_file = 'Pancreas_Immune_cell.h5ad'

adata_concat.obs=adata_concat.obs.drop(labels=['Coarse_Cell_Annotations'],axis=1)
adata_concat.obs['datasets']=adata_concat.obs['batch']
adata_concat.obs=adata_concat.obs.drop(labels=['batch'],axis=1)
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
#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata')
results_file = 'Pancreatic_Immune_cell.h5ad'
adata = sc.read_h5ad(results_file)
adata.obs_names_make_unique()
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

p1 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig('total_counts_genes_scatter2.pdf',format="pdf")
plt.clf()

adata.layers["counts"] = adata.X

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
model_low = models.Model.load(model="Immune_All_Low.pkl")
#duplicated_obs_names = adata.obs_names.duplicated(keep='first')
#adata = adata[~duplicated_obs_names]

adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform

predictions_low = celltypist.annotate(adata_celltypist, model=model_low, majority_voting=True)
predictions_low_adata = predictions_low.to_adata()

adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[adata.obs.index, "majority_voting"]
value_counts = adata.obs["celltypist_cell_label_fine"].value_counts()
adata=reshape_anno(adata)
adata.obs["louvain_1_anno"].value_counts()
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_pancreac_celltypist_cell_label_fine')
adata = adata[adata.obs["louvain_1_anno"] != "Epithelial cells"]
adata = adata[adata.obs["louvain_1_anno"] != "Fibroblasts"]
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)



adata.write(results_file)

def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)
```
