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
from anndata import AnnData
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
sc.set_figure_params(dpi=1000, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
warnings.filterwarnings("ignore")
print(os.getcwd())

#/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma
#########################
#GSE156728
#############
os.chdir('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma')

adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma/UCEC_GSE156728/adata.h5ad')
adata1.obs['patients']=adata1.obs['patient']
adata1.obs['samples']=adata1.obs['patient']
adata1.obs=adata1.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma/UCEC_GSE154763/adata.h5ad')
adata2.obs['patients']=adata2.obs['patient']
adata2.obs['samples']=adata2.obs['patient']
adata2.obs=adata2.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch', 'library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster', 'source', 'tech', 'UMAP1', 'UMAP2'],axis=1)

adata_concat = adata1.concatenate(adata2, batch_categories=['GSE156728','GSE154763']) #合并数据集
sc.tl.pca(adata_concat)
sc.external.pp.bbknn(adata_concat, batch_key='batch')
adata_concat.obs['datasets']=adata_concat.obs.batch
sc.tl.pca(adata_concat, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='EndometrialCarcinoma'
n_comps=50
sc.pl.pca(adata_concat, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')
results_file = "/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/EndometrialCarcinoma_Immune_cell.h5ad"
adata_concat.write(results_file)




adata3=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/cancerscem/UCEC_cancerscem.h5ad')

#var_names = adata1.var_names.union(adata2.var_names) #
#adata1=re_size_adata(adata1,var_names)
#adata2=re_size_adata(adata2,var_names)

#adata = adata1.concatenate(adata2,batch_categories=['GSE156728', 'GSE154763'])
#adata = adata.concatenate(adata3,index_unique=None)
#sc.tl.pca(adata)
#sc.external.pp.bbknn(adata, batch_key='batch')
adata1.obs['datasets']='GSE156728'
#adata.obs=adata.obs.drop(labels=['batch'],axis=1)

#results_file = "/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/EndometrialCarcinoma_integrate_cell.h5ad"
#adata.write(results_file)

results_file = "/share/pub/dengcy/Cancer_Gwas/Runtime1.0/1.data_integrate_progress/integrate_singlecell/EndometrialCarcinoma_Immune_cell.h5ad"
adata1.write(results_file)
```

### 单独处理数据

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


BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma/'
os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file=os.path.join(BaseDirectory, 'adata_pp_immune.h5ad')

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/EndometrialCarcinoma_Immune_cell.h5ad")

adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

sns.displot(adata.obs["total_counts"], bins=100, kde=False)
plt.savefig('total_counts_distplot.pdf',format="pdf")
plt.clf()

# sc.pl.violin(adata, 'total_counts')
sc.pl.violin(adata, "pct_counts_mt")
plt.savefig('pct_counts_mt_violin.pdf',format="pdf")
plt.clf()

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
plt.savefig('total_counts_genes_scatter.pdf',format="pdf")
plt.clf()


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
sc.external.pp.bbknn(adata, batch_key='batch')
sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
color_map='CMRmap'
cancer='EndometrialCarcinoma'
n_comps=50
sc.pl.pca(adata, color='datasets', color_map=color_map, save='_' + str(cancer) + '_' + str(n_comps) + 'comps_PCA')

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

sc.pl.umap(adata,color=["total_counts", "pct_counts_mt"],save='_counts')

import seaborn as sns
import matplotlib.pyplot as plt

sc.tl.louvain(adata, key_added="louvain_res0_25", resolution=0.25)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=1)

adata.write(results_file)
```

### annotation

```python

marker_genes = {"T-cell":["CD3D","CD3E","CD3G"],
                "NKT":["CCL5","NKG7","KLRD1","KLRF1","KLRB1","FGFBP2"],
                "B-cell":["CD79A","CD79B","MS4A1","IGHM","CD19"],
                "Plasma":["IGHG3","IGHG1","MZB1","JCHAIN"],
                "Mast": "TPSAB1",
                "myeloid":["CST3","LYZ","CD68","MS4A6A"],
                "Neutrophil":"S100A8",
                "Fibroblasts":["DCN","COL1A2","COL1A1","COL3A1","ACTA2"],
                "Endothelial":["CLDN5","RAMP2","VWF","CDH5","PECAM1","ENG","FLT1", "CDH1"],
                "EPCAM":["EPCAM","KRT17","KRT5","KRT19","KRT18","KRT15"]
                }
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found

sc.pl.dotplot(adata,groupby="louvain_res1",var_names=marker_genes_in_data,standard_scale="var")
plt.savefig('first_louvain_marker_dotplot.pdf', format="pdf")
plt.clf()


results_file=os.path.join(BaseDirectory, 'adata_pp_immune.h5ad')
adata=sc.read_h5ad(results_file)
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

sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_celltypist_cell_label_fine')

adata.write(results_file)
```