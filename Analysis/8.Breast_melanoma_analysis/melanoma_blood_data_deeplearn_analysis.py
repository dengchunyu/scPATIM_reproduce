import scanpy as sc
import pandas as pd
import seaborn as sns
import os
import numpy as np
import anndata as ad
import torch
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/')
sc.set_figure_params(fontsize=10, dpi=80, dpi_save=300, format='svg')

adata = sc.read_h5ad('melanoma_ICI_blood.h5ad')
adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].astype('category')

adata.obs["first_celltype_annotation"] = adata.obs.louvain_1_anno.map(cl_annotation)

adata.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation'].astype('category')

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="Melanoma ICI pre-treatment in blood",palette=color_annotation,size=6,save='melanoma_blood_umap')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/GSE166181_blood_pre.adata")

adata = adata[(adata.obs['first_celltype_annotation'] == 'Tem/Temra.cyt')]

sc.pl.umap(adata, color='first_celltype_annotation', legend_loc="on data", title="HCCs of Melanoma ICI pre-treatment in blood",palette=color_annotation,size=6,save='melanoma_blood_HCCs_umap')

adata.write("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/GSE166181_blood_pre_HCCs.adata")

adata_tnbc = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_blood_T_pre.adata')

pcc_gene = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/Melanoma/finngen_r7_c3_melanoma_gene_PCC.csv", index_col=0)

# 筛选基因
genes_filter = pcc_gene[(pcc_gene['PCC'] > 0) & (pcc_gene['adj_pvalue'] < 0.01)].index

# 提取500个基因的数据
genes_filter_in_adata = [gene for gene in genes_filter if gene in adata.var_names]

genes_filter_in_adata = [gene for gene in genes_filter_in_adata if gene in adata_tnbc.var_names]

adata = adata[:, genes_filter_in_adata]
adata_tnbc = adata_tnbc[:, genes_filter_in_adata]


# 将分类变量转换为数值
le = LabelEncoder()
adata.obs['response_num'] = le.fit_transform(adata.obs['response'])
adata_tnbc.obs['response_num'] = le.fit_transform(adata_tnbc.obs['response'])
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/')

train_loader, test_loader, feature_importance = MLP_model_test(adata,data_index="GSE166181")
train_loader2, test_loader2, feature_importance2 = MLP_model_test(adata_tnbc,data_index="GSE166181_tnbc_")


### 获得最佳基因
import matplotlib.pyplot as plt
import numpy as np
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE166181/')
hidden_size1 = 256
hidden_size2 = 128
output_size = 2

def get_gene_performance(train_loader, test_loader, feature_importance,data_index):
    sorted_indices = np.argsort(feature_importance)
    performances = []
    best_auc = 0.0
    best_model = None
    def train_and_validate_model(train_loader, test_loader, feature_indices):
        model = DynamicMLP(len(feature_indices), hidden_size1, hidden_size2, output_size)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        criterion = torch.nn.CrossEntropyLoss()
        model.train()
        for inputs, labels in train_loader:
            inputs = inputs[:, feature_indices]
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
        model.eval()
        probs = []
        true_labels = []
        with torch.no_grad():
            for inputs, labels in test_loader:
                inputs = inputs[:, feature_indices]
                outputs = model(inputs)
                prob = torch.softmax(outputs.data, dim=1)[:, 1]  # 假设正类的标签是1
                probs.extend(prob.numpy())
                true_labels.extend(labels.numpy())# 计算ROC曲线
        fpr, tpr, _ = roc_curve(true_labels, probs)
        roc_auc = auc(fpr, tpr)
        return roc_auc, model
    gene_groups = [sorted_indices[:i] for i in range(5, len(sorted_indices) + 1, 5)]
    performances = []
    for group in gene_groups:
        performance, model = train_and_validate_model(train_loader, test_loader, group)
        performances.append(performance)
        if performance > best_auc:
            best_auc = performance
            best_model = model
    plt.clf()
    plt.figure()
    plt.plot(range(5, len(sorted_indices) + 1, 5), performances,color="#EE7214")
    plt.xlabel('Number of Features')
    plt.ylabel('Performance: auc')
    plt.show()
    plt.savefig(data_index+'feature_num_performance.pdf')
    data = {'performance': performances}
    df = pd.DataFrame(data)
    df.to_csv(data_index + 'performances.csv', index=False)

get_gene_performance(train_loader, test_loader,feature_importance,data_index="GSE166181")

get_gene_performance(train_loader2, test_loader2,feature_importance2,data_index="GSE166181_tnbc_")
