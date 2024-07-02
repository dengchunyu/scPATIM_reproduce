import scanpy as sc
import pandas as pd
import scdrs
import numpy as np
import warnings
warnings.filterwarnings("ignore")
adata_file = '/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad'
scDRS_gene_file = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/_gene_PCC.csv'
n_ctrl = 500

adata = sc.read_h5ad(adata_file)
scdrs.preprocess(adata)
if scDRS_gene_file.endswith('.txt'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep='\t')
elif scDRS_gene_file.endswith('.csv'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep=',')

#判断scDRS_gene_file是否有gene，如果没有判断是否有gene_symbol，如果有gene_symbol则将其改为gene，如果没有则报错
if 'GENE' not in scDRS_gene.columns:
    if 'gene_symbol' in scDRS_gene.columns:
        scDRS_gene = scDRS_gene.rename(columns={'gene_symbol':'GENE'})
    else:
        scDRS_gene = scDRS_gene.rename(columns={'Unnamed: 0':'GENE'})
if 'weight_pcc' in scDRS_gene.columns:
    weight_col = 'weight_pcc'

scDRS_gene = scDRS_gene.rename(columns={weight_col:'weight'})
scDRS_gene = scDRS_gene.dropna(subset=['weight'])
max_weight  = scDRS_gene['weight'].replace([np.inf, -np.inf], np.nan).max()
scDRS_gene['weight'] = scDRS_gene['weight'].replace([np.inf, -np.inf],max_weight + 1)

for top_gene_num in range(50, 1050, 50):
    top_gene = scDRS_gene.sort_values(by='weight', ascending=False).iloc[:top_gene_num, :]
    df_score = scdrs.score_cell(data=adata,
                                gene_list=top_gene['GENE'],
                                gene_weight=top_gene['weight'],
                                ctrl_match_key="mean_var",
                                n_ctrl=n_ctrl,
                                weight_opt="vs",
                                return_ctrl_raw_score=False,
                                return_ctrl_norm_score=True,
                                verbose=False)
    scDRS_score_file = f"/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/scDRS_score_top_{top_gene_num}_genes.csv"
    df_score.to_csv(scDRS_score_file,index=False)


import scanpy as sc
import pandas as pd
import scdrs
import numpy as np
import warnings
warnings.filterwarnings("ignore")
adata_file = '/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad'
scDRS_gene_file = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls/_gene_PCC.csv'
n_ctrl = 500

adata = sc.read_h5ad(adata_file)
scdrs.preprocess(adata)
if scDRS_gene_file.endswith('.txt'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep='\t')
elif scDRS_gene_file.endswith('.csv'):
    scDRS_gene = pd.read_csv(scDRS_gene_file, sep=',')

#判断scDRS_gene_file是否有gene，如果没有判断是否有gene_symbol，如果有gene_symbol则将其改为gene，如果没有则报错
if 'GENE' not in scDRS_gene.columns:
    if 'gene_symbol' in scDRS_gene.columns:
        scDRS_gene = scDRS_gene.rename(columns={'gene_symbol':'GENE'})
    else:
        scDRS_gene = scDRS_gene.rename(columns={'Unnamed: 0':'GENE'})
if 'weight_pcc' in scDRS_gene.columns:
    weight_col = 'weight_pcc'

scDRS_gene = scDRS_gene.rename(columns={weight_col:'weight'})
scDRS_gene = scDRS_gene.dropna(subset=['weight'])
max_weight  = scDRS_gene['weight'].replace([np.inf, -np.inf], np.nan).max()
scDRS_gene['weight'] = scDRS_gene['weight'].replace([np.inf, -np.inf],max_weight + 1)

for top_gene_num in range(50, 1050, 50):
    scDRS_gene = scDRS_gene.sort_values(by='weight', ascending=False).iloc[:top_gene_num, :]
    df_score = scdrs.score_cell(data=adata,
                                gene_list=scDRS_gene['GENE'],
                                gene_weight=scDRS_gene['weight'],
                                ctrl_match_key="mean_var",
                                n_ctrl=n_ctrl,
                                weight_opt="vs",
                                return_ctrl_raw_score=False,
                                return_ctrl_norm_score=True,
                                verbose=False)
    scDRS_score_file = f"/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls/scDRS_score_top_{top_gene_num}_genes.csv"
    df_score.to_csv(scDRS_score_file,index=False)

###下面进行可视化
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result")
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")
adata1.obs=adata1.obs.drop(labels=['scPagwas_eqtls_norm_score', 'scPagwas_eqtls_mc_pval', 'scPagwas_norm_score', 'scPagwas_mc_pval', 'scDRS_norm_score', 'scDRS_mc_pval',],axis=1)
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/DLBCL_B_cells.h5ad")
adata2.obs=adata2.obs.drop(labels=['scHII_mc_pval', 'scHII_norm_score', 'scHII_down_mc_pval', 'scHII_down_norm_score', 'scPagwas_down_mc_pval', 'scPagwas_down_norm_score', 'scPagwas_mc_pval', 'scPagwas_norm_score', 'magma_mc_pval', 'magma_norm_score'],axis=1)
auc_results1 = []
auc_results2 = []
for top_gene_num in range(50, 1050, 50):
    scDRS_score_file = f"/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls/scDRS_score_top_{top_gene_num}_genes.csv"
    print(top_gene_num)
    scPatim_score = pd.read_csv(scDRS_score_file)
    scPatim_score.index = adata1.obs.index
    scPatim_score = scPatim_score.iloc[:, 1:3]
    scPatim_score.columns = ['scPatim_' + i for i in scPatim_score.columns]
    if 'scPatim_norm_score' in adata1.obs.columns and 'scPatim_mc_pval' in adata1.obs.columns:
        adata1.obs = adata1.obs.drop(labels=['scPatim_norm_score', 'scPatim_mc_pval'], axis=1)
    adata1.obs = pd.concat([adata1.obs,scPatim_score], axis=1)
    auc=roc_auc_score(adata1.obs['label'],adata1.obs['scPatim_norm_score'])
    auc_results1.append(auc)
    scDRS_score_file = f"/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/scDRS_score_top_{top_gene_num}_genes.csv"
    print(top_gene_num)
    scPatim_score = pd.read_csv(scDRS_score_file)
    scPatim_score.index = adata2.obs.index
    scPatim_score = scPatim_score.iloc[:, 1:3]
    scPatim_score.columns = ['scPatim_' + i for i in scPatim_score.columns]
    if 'scPatim_norm_score' in adata2.obs.columns and 'scPatim_mc_pval' in adata2.obs.columns:
        adata2.obs = adata2.obs.drop(labels=['scPatim_norm_score', 'scPatim_mc_pval'], axis=1)
    adata2.obs = pd.concat([adata2.obs,scPatim_score], axis=1)
    auc=roc_auc_score(adata2.obs['label'],adata2.obs['scPatim_norm_score'])
    auc_results2.append(auc)

import matplotlib.pyplot as plt
data = {'scPATIM_PBMC_auc': auc_results1,
    'scPATIM_DLBLC_auc': auc_results2,
    'top_gene_number': range(50, 1050, 50)}
df = pd.DataFrame(data)
def plot_line(df):
    plt.figure(figsize=(8, 6))
    plt.plot(df['top_gene_number'], df['scPATIM_PBMC_auc'], marker='o', color='#DD5746',label='AUC for PBMC data')
    plt.plot(df['top_gene_number'], df['scPATIM_DLBLC_auc'], marker='s', color='#4793AF',label='AUC for B cells in DLBLC data')
    plt.axvline(x=500, color='gray', linestyle='--',label='500')
    plt.xlabel('Top TRGs numbers are used to calculate the TRS score.')
    plt.ylabel('AUC')
    plt.title('Validation of the top TRGs in two benchmark datasets.')
    plt.legend()
    plt.grid(True)
    plt.savefig("Top_gene_test_for_DLBLC_B.pdf")
    plt.close()

plot_line(df=df)