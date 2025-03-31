import numpy as np
import anndata as ad
import scanpy as sc
import pandas as pd
import seaborn as sns
import pyarrow as pa
import os
import pyarrow.feather as feather
from scipy.sparse import coo_matrix
import scipy.sparse
import csv

# 设置随机种子
np.random.seed(42)

### 1. 加载数据
# 加载 GWAS 数据
gwas_data = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt", sep="\t")

# 加载单细胞数据
sc_data = ad.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")

### 2. 添加噪声的函数
def add_noise_to_gwas(gwas_df, noise_level):
    """
    给GWAS数据的beta列添加噪声。
    :param gwas_df: GWAS数据的DataFrame
    :param noise_level: 噪声比例 (k)
    :return: 带噪声的GWAS数据
    """
    beta_col = gwas_df["beta"]
    noise = np.random.normal(0, noise_level * np.std(beta_col), size=beta_col.shape)
    gwas_df["beta_noisy"] = beta_col + noise
    return gwas_df

def add_noise_to_sc_data(sc_adata, noise_level):
    """
    给单细胞数据的基因表达矩阵添加噪声。
    :param sc_adata: 单细胞数据的AnnData对象
    :param noise_level: 噪声比例 (k)
    :return: 带噪声的单细胞数据
    """
    noisy_X = sc_adata.X + np.random.normal(0, noise_level * np.mean(sc_adata.X), sc_adata.X.shape)
    noisy_adata = sc_adata.copy()
    noisy_adata.X = np.clip(noisy_X, 0, None)
    return noisy_adata

### 4. 噪声水平测试
noise_levels = [0.1, 0.2, 0.3, 0.4, 0.5]  # 噪声比例
BaseDirectory="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/"
for noise_level in noise_levels:
    noisy_gwas = add_noise_to_gwas(gwas_data.copy(), noise_level)
    noisy_sc_data = add_noise_to_sc_data(sc_data, noise_level)
    noisy_gwas.to_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/gwas_noise" + str(noise_level) + ".csv")
    noisy_sc_data.write("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scdata_noise" + str(noise_level) + ".h5ad")
    if scipy.sparse.issparse(noisy_sc_data.X):
        dense_matrix = noisy_sc_data.X.toarray()
    else:
        dense_matrix = noisy_sc_data.X
    dense_matrix = np.transpose(dense_matrix)
    df = pd.DataFrame(dense_matrix)
    obs=noisy_sc_data.obs
    df.index = noisy_sc_data.var_names
    df.columns=obs.index
    out_file=os.path.join(BaseDirectory,'matrix_'+str(noise_level) +'.feather')
    feather.write_feather(df, out_file)
    obs=noisy_sc_data.obs
    obs_file=os.path.join(BaseDirectory,'obs'+str(noise_level) + '.csv')
    obs.to_csv(obs_file, index=False, header=True)
    var=noisy_sc_data.var_names
    var_file=os.path.join(BaseDirectory,'var_names'+str(noise_level) + '.csv')
    var=pd.DataFrame(var)
    var.to_csv(var_file, index=False, header=True)



#!/bin/bash
#SBATCH -e noise.err
#SBATCH -o noise.out
#SBATCH -J noise
#SBATCH -w in008
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
python /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/noise_data.py