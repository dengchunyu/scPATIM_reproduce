import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/GSE120575.h5ad")
adata.obs.response.value_counts()
#Non-responder    10190
#Responder         5110
adata.obs.therapy.value_counts()
#anti-PD1          11182
#anti-CTLA4+PD1     3601
#anti-CTLA4          517
adata.obs.post_pre.value_counts()
#Post    9372
#Pre     5928
adata.obs.celltype.value_counts()
#celltype
#Exhausted CD8 T cells         6128
#Regulatory T cells            3373
#Cytotoxicity (Lymphocytes)    1818
#B cells                       1682
#Monocytes/Macrophages         1394
#Memory T cells                 624
#Dendritic cells                281
#Name: count, dtype: int64
import celltypist
from celltypist import models
adata = adata[adata.obs['post_pre'] == 'Pre']
adata.layers["counts"] = adata.X
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)

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
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].astype('category')
adata.obs["first_celltype_annotation"] = adata.obs.louvain_1_anno.map(cl_annotation)
adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].astype('category')
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/6.therapy_analysis")

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="melanoma ICI pre-treatment",palette=color_annotation,size=5,add_outline=True, save='melanoma_singlecell_umap')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/GSE120575_tissue_pre.adata")
adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/GSE120575_tissue_pre.adata")
# 提取特定细胞类型
adata_hccs = adata[(adata.obs['first_celltype_annotation'] == 'Tem/Trm.cyt') | (adata.obs['first_celltype_annotation'] == 'Tem/Temra.cyt')]

sc.pl.umap(adata_hccs, color='first_celltype_annotation', legend_loc="on data", title="melanoma ICI pre-treatment in HCCs",palette=color_annotation,size=5,add_outline=True, save='melanoma_HCCs_umap')

celltype_counts = specific_pre_cells.obs['celltype'].value_counts()
response_counts = specific_pre_cells.obs['response'].value_counts()

# 对每个celltype和response类别进行计数
celltype_response_counts = specific_pre_cells.obs.groupby(['celltype', 'response']).size()

adata_hccs.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/GSE120575_tissue_pre_HCCs.h5ad")

