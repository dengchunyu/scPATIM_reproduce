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
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#
cancers=['Pancreatic','BreastCancer','HeadandNeck','EsophagealCancer','ThyroidCancer','BladderCancer','ColorectalCancer','LiverCancer','LungCancer','Melanoma','OvarianCancer','EndometrialCarcinoma','KidneyCancer','ProstateCancer','GastricCancer']
os.chdir(BaseDirectory)
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='svg')
##给不同细胞类型搞一个缩写形式，描写方便
#查看共有多少细胞类型
# 创建一个空的数据框来存储所有癌症的统计结果
combined_stat_df = pd.DataFrame()
celltype_list=[]
# 迭代处理每个癌症数据
#    samples = len(obs['samples'].unique())
#    patients = len(obs['patients'].unique())
#    Ndatasets = len(obs['datasets'].unique())
for i in cancers:
    print(i)
    results_file = os.path.join(BaseDirectory, i + '_Immune_cell.h5ad')
    adata = sc.read_h5ad(results_file)
    celltype_list.append(adata.obs['louvain_1_anno'].unique())
    obs = adata.obs
    cellnumber = len(obs['louvain_1_anno'])
    Ncelltype = len(obs['louvain_1_anno'].unique())
    stat_df = pd.DataFrame({'cancer': [i], 'cellnumber': [cellnumber], 'Ncelltype': [Ncelltype]})
    combined_stat_df = pd.concat([combined_stat_df, stat_df], ignore_index=True)

# 保存整合的统计结果到一个文件
combined_stat_df.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/combined_stat_file.csv', index=False, sep=',')
merged_list = [cell for sublist in celltype_list for cell in sublist]
unique_cells = set(merged_list)