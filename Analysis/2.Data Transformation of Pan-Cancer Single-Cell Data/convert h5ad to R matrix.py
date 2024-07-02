import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import pyarrow as pa
import os
import pyarrow.feather as feather
from scipy.sparse import coo_matrix
import scipy.sparse
import csv
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#
cancers=['ProstateCancer','ColorectalCancer','LiverCancer','LungCancer','HeadandNeck','Melanoma','OvarianCancer','Pancreatic','GastricCancer','BreastCancer','EndometrialCarcinoma','KidneyCancer','ThyroidCancer','EsophagealCancer']
#i=,
for i in cancers:
    print(i)
    results_file='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'+i+'_Immune_cell.h5ad'
    adata=sc.read_h5ad(results_file)
    if scipy.sparse.issparse(adata.X):
        dense_matrix = adata.X.toarray()
    else:
        dense_matrix = adata.X
    dense_matrix = np.transpose(dense_matrix)
    df = pd.DataFrame(dense_matrix)
    obs=adata.obs
    df.index = adata.var_names
    df.columns=obs.index
    out_file=os.path.join(BaseDirectory,i+'_matrix.feather')
    feather.write_feather(df, out_file)
    obs=adata.obs
    obs_file=os.path.join(BaseDirectory,i+'_obs.csv')
    obs.to_csv(obs_file, index=False, header=True)
    var=adata.var_names
    var_file=os.path.join(BaseDirectory,i+'_var_names.csv')
    var=pd.DataFrame(var)
    var.to_csv(var_file, index=False, header=True)

for i in cancers:
    print(i)
    results_file='/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'+i+'_Immune_cell.h5ad'
    adata=sc.read_h5ad(results_file)
    obs=adata.obs
    obs_file=os.path.join(BaseDirectory,i+'_obs.csv')
    if 'merge_celltype_annotation' in obs.columns:
        obs.to_csv(obs_file, index=False)
        print("Merge列存在，并已成功输出到CSV文件。")
    else:
        print("Merge列不存在。")