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
from scipy.stats import norm

BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#
cancers=['HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")


for i in cancers:
    print(i)
    results_file=os.path.join(BaseDirectory,i+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    obs=adata.obs
    ct_value=merge_celltype_p(obs['scDRS_pval'], obs['first_celltype_annotation'])
    ct_value.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/'+i+'/first_celltype_pvalue.csv')
    ct_value=merge_celltype_p(obs['scDRS_pval'], obs['merge_celltype_annotation'])
    ct_value.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/'+i+'/Merged_celltype_pvalue.csv')
    adata.obs['Significant_cells'] = adata.obs['scDRS_pval'].apply(lambda x: 1 if x < 0.05 else 0)
    grouped = obs.groupby('merge_celltype_annotation')['Significant_cells']
    grouped_percentage = (grouped.sum() / grouped.count()) * 100
    average_scores = obs.groupby('merge_celltype_annotation')['scDRS_norm_score'].mean()
    ct_value['sig_percentage']=grouped_percentage[ct_value['celltype']].values
    ct_value['average_trs']=average_scores[ct_value['celltype']].values
    output_file = f'{i}_merge_average_pscores.csv'
    ct_value.to_csv(output_file, index=False)

def p_merge(pvalues):
    zvalues = -np.sqrt(2) * norm.ppf(pvalues / 2)
    ztotal = np.mean(zvalues)
    p_total = norm.cdf(-abs(ztotal))
    return p_total

def merge_celltype_p(single_p, celltype):
    celltype_p = pd.DataFrame({'celltype': celltype, 'pvalue': single_p})
    celltype_p = celltype_p.groupby('celltype')['pvalue'].agg(p_merge).reset_index()
    return celltype_p
