# bladdercancer

/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/BladderCancer_data/GSE145281

## 预处理数据

```python
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/BladderCancer_data/GSE145281")

umi_matrix_file = "GSE145281_All_responder_nonresponder_raw_matrix.txt"
adata = sc.read_text(umi_matrix_file)
adata=adata.T

sample_response_df = pd.DataFrame(index=adata.obs_names, columns=['sample', 'response'])
# 根据adata.obs_names设置'sample'和'response'
sample_response_df['sample'] = sample_response_df.index.str.split('_').str[0]
sample_response_df['response'] = 'response'
sample_response_df.loc[sample_response_df['sample'].str.startswith('NR'), 'response'] = 'non-response'
# 将'sample'和'response'合并到adata.obs中
adata.obs = pd.concat([adata.obs, sample_response_df], axis=1)

adata.write("GSE145281_ici.h5ad")

def progress_f(adata):
    adata.layers["counts"] = adata.X
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata, key_added="louvain_res1", resolution=0.5)
    return adata
adata=progress_f(adata)
adata.write("GSE145281_ici.h5ad")


def reshape_anno(adata):
    sub_df = adata.obs[['louvain_res1', 'celltypist_cell_label_fine']]
    pivot_table = pd.pivot_table(sub_df, index='louvain_res1', columns='celltypist_cell_label_fine', aggfunc='size', fill_value=0)
    celltype_names = pivot_table.idxmax(axis=1)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_res1'].map(celltype_names)
    adata.obs['louvain_1_anno'] = adata.obs['louvain_1_anno'].cat.remove_unused_categories()
    return(adata)
import celltypist
from celltypist import models
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
adata = adata[(adata.obs['louvain_1_anno'] == 'Naive B cells')]
adata.write("GSE145281_ici.h5ad")
```

## 深度学习训练遗传基因获得免疫应答marker

```python
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/BladderCancer_data/GSE145281")

adata = sc.read_h5ad("GSE145281_ici.h5ad")
pcc_gene = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/BladderCancer/finngen_r7_c3_bladder_exallc_gene_PCC.csv", index_col=0)
# 筛选基因
genes_filter = pcc_gene[(pcc_gene['PCC'] > 0) & (pcc_gene['adj_pvalue'] < 0.01)].index
#基因数量
#len(genes_filter)
#136

import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/BladderCancer_data/GSE145281")

# 提取基因的数据
genes_filter_in_adata = [gene for gene in genes_filter if gene in adata.var_names]
adata = adata[:, genes_filter_in_adata]

# 将分类变量转换为数值
le = LabelEncoder()
adata.obs['response_num'] = le.fit_transform(adata.obs['response'])

#运行训练模型
train_loader, test_loader, feature_importance = MLP_model_test(adata,data_index="GSE145281")

# 加载模型
model = DynamicMLP(input_size, hidden_size1, hidden_size2, output_size)  # 重新创建一个模型实例
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)  # 重新创建一个优化器实例

checkpoint = torch.load('model_checkpoint.pth')
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
```

### 获得最佳基因

```python

import matplotlib.pyplot as plt
import numpy as np

from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score

#你的模型是一个二分类问题，所以你可以使用torch.nn.CrossEntropyLoss作为你的损失函数。你可以在初始化模型和优化器的部分添加如下代码来定义你的损失函数
input_size = len(genes_filter_in_adata)
hidden_size1 = 60
hidden_size1 = 30
output_size = 2
sorted_indices = np.argsort(feature_importance)
def train_and_validate_model(train_loader, test_loader, feature_importance):
    sorted_indices = np.argsort(feature_importance)
    performances = []
    best_auc = 0.0
    best_model = None
    model = DynamicMLP(len(sorted_indices), hidden_size1, hidden_size2, output_size)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = torch.nn.CrossEntropyLoss()
    model.train()
    for inputs, labels in train_loader:
        inputs = inputs[:, sorted_indices]
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
            inputs = inputs[:, sorted_indices]
            outputs = model(inputs)
            prob = torch.softmax(outputs.data, dim=1)[:, 1]  # 假设正类的标签是1
            probs.extend(prob.numpy())
            true_labels.extend(labels.numpy())# 计算ROC曲线
    fpr, tpr, _ = roc_curve(true_labels, probs)
    roc_auc = auc(fpr, tpr)
    return roc_auc, model

# 遍历每个特征数量
gene_groups = [sorted_indices[:i] for i in range(5, len(sorted_indices) + 1, 5)]
performances = []
for group in gene_groups:
    performance, model = train_and_validate_model(train_loader, test_loader, group)
    performances.append(performance)
    if performance > best_auc:
        best_auc = performance
        best_model = model

# 绘制特征数量与性能的关系图
def plot_ferformance(sorted_indices,performances):
    plt.clf()
    plt.figure()
    plt.plot(range(5, len(sorted_indices) + 1, 5), performances,color="#EE7214")
    plt.xlabel('Number of Features')
    plt.ylabel('Performance: auc')
    plt.show()
    plt.savefig('feature_num_performance.pdf')
    data = {'performance': performances}
    df = pd.DataFrame(data)
    df.to_csv('performances.csv', index=False)

plot_ferformance(sorted_indices,performances)
```
函数
```python
# 定义Dataset
class GeneDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)  # 转换为稠密矩阵
        self.y = torch.tensor(y, dtype=torch.long)
    def __len__(self):
        return self.X.shape[0]  # 使用稠密矩阵的行数作为长度
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]


# 定义MLP模型
class DynamicMLP(nn.Module):
    def __init__(self, input_size, hidden_size1, hidden_size2, output_size):
        super(DynamicMLP, self).__init__()
        self.layers = nn.Sequential(
            nn.Linear(input_size, hidden_size1),
            nn.ReLU(),
            nn.Linear(hidden_size1, hidden_size2),
            nn.ReLU(),
            nn.Linear(hidden_size2, output_size)
        )
    def forward(self, x):
        x = self.layers(x)
        return x

def MLP_model_test(adata,data_index):
    X_train, X_test, y_train, y_test = train_test_split(adata.X, adata.obs['response_num'], test_size=0.2, random_state=42)
    train_data = GeneDataset(X_train, y_train)
    test_data = GeneDataset(X_test, y_test)
    input_size = len(genes_filter_in_adata)
    hidden_size1 = 60
    hidden_size1 = 30
    output_size = 2
    train_data = GeneDataset(X_train, y_train)
    test_data = GeneDataset(X_test, y_test)
    train_loader = DataLoader(train_data, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_data, batch_size=32, shuffle=False)
    model = DynamicMLP(input_size, hidden_size1, hidden_size2, output_size)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(100):
        for i, (inputs, labels) in enumerate(train_loader):
            outputs = model(inputs)
            loss = nn.CrossEntropyLoss()(outputs, labels)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
    model.eval()
    probs = []
    true_labels = []
    with torch.no_grad():
        for inputs, labels in test_loader:
            outputs = model(inputs)
            prob = torch.softmax(outputs.data, dim=1)[:, 1]  # 假设正类的标签是1
            probs.extend(prob.numpy())
            true_labels.extend(labels.numpy())
    fpr, tpr, _ = roc_curve(true_labels, probs)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='#A25772', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='#7C93C3', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.savefig(data_index + '_roc_curve.pdf')
    torch.save({'model_state_dict': model.state_dict(),'optimizer_state_dict': optimizer.state_dict()}, 'model_checkpoint.pth')
    grads = torch.zeros((input_size,))
    for inputs, labels in train_loader:
        inputs.requires_grad = True  # 设置requires_grad=True以计算梯度# 执行前向传播
        outputs = model(inputs)# 选择一个类别计算梯度
        output_idx = outputs.argmax(dim=1)
        output_max = outputs[torch.arange(outputs.shape[0]), output_idx]# 执行反向传播，计算梯度
        output_max.backward(torch.ones_like(output_max))# 累积每个批次的梯度
        grads += inputs.grad.abs().sum(dim=0).detach()# 清空梯度，为下一轮计算做准备
        model.zero_grad()
    feature_importance = grads / len(train_loader.dataset)
    print(feature_importance.argsort(descending=True)[:10])
    important_genes_indices = feature_importance.argsort(descending=True)
    important_genes = [genes_filter_in_adata[i] for i in important_genes_indices]
    data = {'important_genes': important_genes, 'feature_importance': feature_importance[important_genes_indices]}
    df = pd.DataFrame(data)
    df.to_csv("feature_importance.csv")
    plt.clf()
    top_10_genes = df.head(10)
    top_10_genes = top_10_genes.sort_values(by='feature_importance', ascending=True)
    plt.clf()
    plt.figure(figsize=(4, 5))
    plt.barh(top_10_genes['important_genes'], top_10_genes['feature_importance'],color="#DBCC95")
    plt.xlabel('Feature Importance')
    plt.ylabel('Important Genes')
    plt.title('Top 10 Important Genes and Their Feature Importance')
    plt.tight_layout()
    plt.savefig('top_10_genes_importance.pdf')
    return train_loader, test_loader, feature_importance

```
