import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246")

adata = sc.read_mtx("./scRNAseq/matrix.mtx.gz")
adata=adata.transpose()
#读取GSE169246_TNBC_RNA.barcode.tsv.gz作为adata的行名
adata.obs_names = pd.read_csv("./scRNAseq/barcodes.tsv.gz", header=None, sep="\t")[0]
#读取features.tsv.gz作为adata的列名
adata.var_names = pd.read_csv("./scRNAseq/features.tsv.gz", header=None, sep="\t")[0]

#读取obs数据
obs = pd.read_csv("./GSE169246_obs.csv", sep=",", header=0, index_col=False)
obs.index = obs["Cell barcode"]
#查看obs.index和adata.obs_names是否一致
obs.index.isin(adata.obs_names).sum()
if obs.index.isin(adata.obs_names).sum() == obs.shape[0]:
      print("obs.index is the same with adata.obs_names")
      adata = adata[obs.index, :]
      adata.obs = obs
      adata.write("./GSE169246.h5ad")
else:
      print("obs.index is not the same with adata.obs_names")
adata= sc.read_h5ad("./GSE169246.h5ad")
adata_blood = adata[adata.obs['Tissue'] == 'blood']
adata_blood = adata_blood[adata_blood.obs['Group'] == 'Pre-treatment']

adata_tissue = adata[adata.obs['Tissue'] != 'blood']
adata_tissue = adata_tissue[adata_tissue.obs['Tissue'] != 'brain']
adata_tissue = adata_tissue[adata_tissue.obs['Tissue'] != 'lung']
adata_tissue = adata_tissue[adata_tissue.obs['Group'] == 'Pre-treatment']
adata_blood.write("GSE169246_blood_pre.adata")
adata_tissue.write("GSE169246_tissue_pre.adata")

###处理并且注释
import celltypist
from celltypist import models

adata = sc.read_h5ad("GSE169246_blood_pre.adata")
adata = sc.read_h5ad("GSE169246_tissue_pre.adata")


adata.layers["counts"] = adata.X
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

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
adata.write("GSE169246_blood_pre.adata")
adata.write("GSE169246_tissue_pre.adata")

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_tissue_pre.adata")
cl_annotation = {'DC':'DC', 'B cells':'Bcells' , 'Regulatory T cells':'Treg', 'Germinal center B cells':'Bgc', 'Tem/Temra cytotoxic T cells':'Tem/Temra.cyt', 'CD16- NK cells':'NK.CD16neg', 'Proliferative germinal center B cells':'Bpgc', 'Tcm/Naive cytotoxic T cells':'Tcm/N.cyt', 'DC2':'DC2', 'pDC':'pDC', 'Alveolar macrophages':'Mac', 'Non-classical monocytes':'Mono.nc', 'Tcm/Naive helper T cells':'Tcm/N.h', 'Tem/Trm cytotoxic T cells':'Tem/Trm.cyt', 'NK cells':'NK', 'Trm cytotoxic T cells':'Trm.cyt', 'Type 17 helper T cells':'Th17', 'Tem/Effector helper T cells':'Tem/Eff.h', 'CD16+ NK cells':'NK.CD16pos', 'Migratory DCs':'DC.mig', 'Naive B cells':'Bnaive', 'Mast cells':'Mast', 'Macrophages':'Mac','Classical monocytes':'Mono.c', 'ILC':'ILC', 'Plasma cells':'Plasma', 'Memory B cells':'Bmem', 'Intermediate macrophages':'Mac.inter', 'gamma-delta T cells':'Tgd','MAIT cells': 'MAIT','Erythrophagocytic macrophages':'Ery.Mac','HSC/MPP':'HSC/MPP'}

color_annotation = {'DC':'#1f77b4', 
  'Bcells' :'#ff7f0e', 
  'Treg':'#279e68', 
  'Bgc':'#023fa5', 
  'Tem/Temra.cyt':"#36600E", 
  'NK.CD16neg':"#9569AB", 
  'Bpgc':'#7d87b9', 
  'Tcm/N.cyt':'#bec1d4', 
  'DC2':"#CA8C74", 
  'pDC':'#bb7784', 
  'Mono.nc':'#8c6d31', 
  'Tcm/N.h':'#ad494a', 
  'Tem/Trm.cyt':'#d62728', 
  'NK': '#f7b6d2', 
  'Trm.cyt':'#dbdb8d', 
  'Th17':'#c49c94', 
  'Tem/Eff.h':'#c5b0d5', 
  'NK.CD16pos':'#ff9896', 
  'DC.mig':'#98df8a', 
  'Bnaive':'#ffbb78', 
  'Mast':'#aec7e8', 
  'Mac':'#17becf',
  'Mono.c':'#b5bd61', 
  'Plasma':'#e377c2', 
  'Bmem':"#E77A77",
  'Mac.inter':'#aa40fc', 
  'Tgd':'#9edae5',
  'MAIT':'#6A8473',
  'Ery.Mac':'#FFD6A5',
  'HSC/MPP':'#4a6fe3'}
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/6.therapy_analysis")
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='svg')

adata = adata[(adata.obs['Treatment'] == 'Anti-PD-L1+Chemo')]
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].astype('category')
adata.obs["first_celltype_annotation"] = adata.obs.louvain_1_anno.map(cl_annotation)
adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].astype('category')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_tissue_pre.adata")

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="BreastCancer ICI pre-treatment",palette=color_annotation,size=2,add_outline=True,save='BreastCancer_singlecell_umap')


adata.obs.first_celltype_annotation.value_counts()
adata = adata[(adata.obs['first_celltype_annotation'] == 'Tem/Trm.cyt')]

#P019,P010,P012,P007
response_patients = ['P019', 'P010', 'P012', 'P007']
adata.obs['response'] = 'non-response'
adata.obs.loc[adata.obs['Patient'].isin(response_patients), 'response'] = 'response'
adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_Tissue_T_pre.adata")
adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_Tissue_T_pre.adata")
sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="BreastCancer HCCs for ICI pre-treatment", palette=color_annotation,add_outline=True,save='BreastCancer_HCCs_umap')


###可视化blood细胞
adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_blood_pre.adata")
adata = adata[(adata.obs['Treatment'] == 'Anti-PD-L1+Chemo')]
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].astype('category')

adata.obs["first_celltype_annotation"] = adata.obs.louvain_1_anno.map(cl_annotation)

adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].astype('category')

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="BreastCancer ICI pre-treatment in blood",palette=color_annotation,size=2,save='BreastCancer_blood_umap')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_blood_pre.adata")


adata = adata[(adata.obs['first_celltype_annotation'] == 'Tem/Temra.cyt') | (adata.obs['first_celltype_annotation'] == 'NK.CD16pos')]
adata.obs['response'] = 'non-response'
adata.obs.loc[adata.obs['Patient'].isin(response_patients), 'response'] = 'response'

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="BreastCancer ICI pre-treatment in blood HCCs",palette=color_annotation,size=2,save='BreastCancer_blood_HCCs_umap')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_blood_T_pre.adata")
