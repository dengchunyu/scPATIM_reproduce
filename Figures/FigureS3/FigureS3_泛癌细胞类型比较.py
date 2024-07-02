import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import numpy as np
import warnings
import gc
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/')
cancer_types = ['HeadandNeck', 'BreastCancer', 'ProstateCancer', 'Melanoma', 'LungCancer',
                'BladderCancer', 'KidneyCancer', 'ThyroidCancer', 'OvarianCancer', 'Pancreatic',
                'GastricCancer', 'LiverCancer', 'ColorectalCancer', 'EsophagealCancer', 'EndometrialCarcinoma']

cancer_data = {}

for cancer_type in cancer_types:
    file_path = f'/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/{cancer_type}_Immune_cell.h5ad'
    adata = sc.read(file_path)
    cancer_data[cancer_type] = adata

all_cancer_data = pd.concat([adata.obs['merge_celltype_annotation'].rename(cancer_type) for cancer_type, adata in cancer_data.items()], axis=1)

# 添加每个细胞的癌症信息
cancer_info = pd.DataFrame({'CancerType': all_cancer_data.columns})

# 将合并的细胞类型信息和癌症信息连接起来
merged_data = pd.melt(all_cancer_data, var_name='CancerType', value_name='CellType')
merged_data = pd.DataFrame(merged_data)  # 替换为你的实际数据

# 删除包含 NaN、Tgd 和 MAIT 的行
merged_data = merged_data.dropna(subset=['CellType'])
merged_data = merged_data[~merged_data['CellType'].isin(['Tgd', 'MAIT'])]

# 计算每个细胞类型在每个癌症中的比例
celltype_counts = pd.crosstab(merged_data['CancerType'], merged_data['CellType'], normalize='index')

color_annotation = {'DC':'#1f77b4', 
  'B/Plasma' :'#ff7f0e', 
  'Treg':'#279e68', 
  'Th':'#c5b0d5',
  'T.cyt':"#36600E", 
  'NK':"#9569AB", 
  'Mac/Mono':'#8c6d31', 
  'Mast':'#aec7e8',
  'Plasma':'#e377c2',
}
cell_types = list(celltype_counts.columns)

# 根据你的数据创建颜色列表
colors = [color_annotation[cell] for cell in cell_types]
import matplotlib.pyplot as plt
# 绘制堆叠柱状图
plt.figure(figsize=(12, 8))
celltype_counts.plot(kind='bar', stacked=True, color=colors)

# 添加标签和标题
plt.ylabel('Percentage of Celltypes')
plt.title('Cell Type Distribution Across All Cancer Types')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(False)
# 显示图形
plt.tight_layout()
plt.show()

# 保存到pdf文件
plt.savefig('percentage_by_cancer_type.pdf')
