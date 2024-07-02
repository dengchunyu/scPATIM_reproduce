import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import mutual_info_score

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.2othergwas_HIIvalid/NM_Healthy_pbmc.h5ad")
trs_df=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls/_singlecell_scPagwas_score_pvalue.Result.csv")
trs_df.index = trs_df.index.astype(str)

# 选择 adata 中与 trs_df.index 对应的标签的子集
adata_sub = adata[adata.obs_names.isin(trs_df.index)].copy()

scPagwas_gPAS_score = trs_df['scPagwas.gPAS.score'].values
data = np.maximum(scPagwas_gPAS_score, 0)  # 保留右峰数据
labels = adata_sub.obs['label'].values
la = pd.Series(labels)
la.value_counts()
#0    88876
#1     8159
#Name: count, dtype: int64
results = {}
for percentile in range(90, 100):
    threshold = np.percentile(data, percentile)
    selected_data = data[data >= threshold]
    selected_labels = labels[data >= threshold]
    if len(selected_data) < 2:
        continue
    clf = DecisionTreeClassifier(max_depth=1)
    clf.fit(selected_data.reshape(-1, 1), selected_labels)
    predicted_labels = clf.predict(selected_data.reshape(-1, 1))
    info_gain = mutual_info_score(selected_labels, predicted_labels)
    results[percentile] = info_gain
best_percentile = max(results, key=results.get)
best_info_gain = results[best_percentile]

print(f'最优分位数为: {best_percentile}th percentile, 信息增益为: {best_info_gain}')
###结果不好

trs_df=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_B_scPATIM/_singlecell_scPagwas_score_pvalue.Result.csv")