import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import numpy as np
import anndata as ad
import matplotlib.pyplot as pt
from matplotlib.pyplot import rc_context,plot,savefig
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/')

matrix_file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/GSE166181_Raw_UMI_CountMatrix_unnormalized_.tsv"
adata = pd.read_csv(matrix_file,sep="\t",index_col=0)
n_obs=len(adata.columns)
n_vars = len(adata.index)
obs = pd.DataFrame()
obs.index=adata.columns
var_names=adata.index
var = pd.DataFrame(index=var_names)
adata=adata.T
adata=ad.AnnData(adata,obs=obs,var=var,dtype='int32')
adata.var_names_make_unique()

anno_data = pd.read_table('GSE166181_Metadata.tsv',index_col=0)
# Append metadata
print('Appending metadata')
adata.obs = anno_data.loc[adata.obs_names]
adata.write("GSE166181.adata")

adata_pre = adata[(adata.obs['time'] == 'NR T0') | (adata.obs['time'] == 'R T0')]
adata_pre.write("GSE166181_pre.adata")

adata = sc.read_h5ad("GSE166181_pre.adata")
sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

results_file='melanoma_ICI_blood.h5ad'
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
from scipy.stats import median_abs_deviation

adata.obs["outlier"] = (is_outlier(adata, "log1p_total_counts", 5) | is_outlier(adata, "log1p_n_genes_by_counts", 5) | is_outlier(adata, "pct_counts_in_top_20_genes", 5))

adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (adata.obs["pct_counts_mt"] > 8)
adata.obs.outlier.value_counts()

adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

adata.layers["counts"] = adata.X

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)

import celltypist
from celltypist import models
model_low = models.Model.load(model="Immune_All_Low.pkl")
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["counts"]  # set adata.X to raw counts
adata_celltypist.X = adata_celltypist.X.astype(float)
sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=float(10000))

sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
#adata_celltypist.X = adata_celltypist.X.toarray()

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
sc.pl.umap(adata,color=["louvain_res1","louvain_1_anno"],save='_melanoma_ICI_blood_umap')
adata.write(results_file)