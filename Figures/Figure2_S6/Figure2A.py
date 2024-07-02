import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import pyarrow as pa
import os
import scipy.sparse
import csv
import matplotlib.pyplot as plt
import matplotlib.colors as colors
##可视化所有癌症类型中显著的细胞类型，不显著的细胞类型显示为灰色
cancers=['LungCancer','BreastCancer', 'Melanoma','GastricCancer','EndometrialCarcinoma','ThyroidCancer','Pancreatic','EsophagealCancer','BladderCancer','ColorectalCancer','LiverCancer','KidneyCancer','OvarianCancer','ProstateCancer','HeadandNeck']

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/sig_celltype_umap')
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='svg')

for i in range(0,15):
    print(i)
    cancer=cancers[i]
    results_file=os.path.join(BaseDirectory,cancer+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    ct_value=pd.read_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/'+cancer+'/first_celltype_pvalue.csv')
    sig = ct_value[ct_value['pvalue'] < 0.05]
    sig_c=sig['celltype'].tolist()
    adata.obs['sig_celltypes'] = adata.obs['first_celltype_annotation'].apply(lambda x: x if x in sig_c else 'Non-significant')
    adata.obs['sig_celltypes'] = adata.obs['sig_celltypes'].astype('category')
    sc.pl.umap(adata, color='sig_celltypes', title=cancer,palette=color_annotation,save=cancer + '_sig_celltype')


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
  'HSC/MPP':'#4a6fe3',
  'Non-significant':'#A9A9A9'}
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='svg')
color_map = color_annotation
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
colors_list = ['#86ABA1', '#DFF3E3', '#F1AE89', '#D68060']
custom_cmap1 = colors.ListedColormap(colors_list)
colors_list = ['#EBE4D1', '#E55604']
custom_cmap2 = colors.ListedColormap(colors_list)
#cancers=['EndometrialCarcinoma','GastricCancer']
#,
cancer_list = ['HeadandNeck','BreastCancer','ThyroidCancer','LungCancer','Melanoma','GastricCancer','Pancreatic','BladderCancer','EsophagealCancer','ColorectalCancer','LiverCancer','OvarianCancer','KidneyCancer','ProstateCancer','EndometrialCarcinoma']

fig, axs = plt.subplots(nrows=15, figsize=(5,80), subplot_kw={'sharex': None, 'sharey': None })

for i in range(0,15):
    print(i)
    cancer=cancer_list[i]
    results_file=os.path.join(BaseDirectory,cancer+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title=cancer, add_outline=True, palette=color_annotation,ax=axs[i])
plt.savefig('cancer_celltype_conbineplot.svg')
