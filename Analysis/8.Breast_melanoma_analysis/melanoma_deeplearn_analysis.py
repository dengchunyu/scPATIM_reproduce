import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

#读取基因
pcc_gene = pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/Melanoma/finngen_r7_c3_melanoma_gene_PCC.csv", index_col=0)

# 筛选基因
genes_filter = pcc_gene[(pcc_gene['PCC'] > 0) & (pcc_gene['adj_pvalue'] < 0.01)].index
#基因数量
#len(genes_filter)
#508
adata=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/GSE120575_tissue_pre_HCCs.h5ad")
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

# 提基因的数据
genes_filter_in_adata = [gene for gene in genes_filter if gene in adata.var_names]

adata_brest=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/TNBC_data/GSE169246/GSE169246_Tissue_T_pre.adata")
genes_filter_in_adata = [gene for gene in genes_filter_in_adata if gene in adata_brest.var_names]

adata = adata[:, genes_filter_in_adata]
adata_brest = adata_brest[:, genes_filter_in_adata]

# 将分类变量转换为数值
le = LabelEncoder()
adata.obs['response_num'] = le.fit_transform(adata.obs['response'])
adata_brest.obs['response_num'] = le.fit_transform(adata_brest.obs['response'])

# 划分训练集和测试集
X_train1, X_test1, y_train1, y_test1 = train_test_split(adata.X, adata.obs['response_num'], test_size=0.2, random_state=42)

X_train2, X_test2, y_train2, y_test2 = train_test_split(adata_brest.X, adata_brest.obs['response_num'], test_size=0.2, random_state=42)

# 定义Dataset
class GeneDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)
    def __len__(self):
        return len(self.y)
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]
class GeneDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X.toarray(), dtype=torch.float32)  # 转换为稠密矩阵
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

# 转换数据为PyTorch数据集
train_data1 = GeneDataset(X_train1, y_train1)
test_data1 = GeneDataset(X_test1, y_test1)

train_data2 = GeneDataset2(X_train2, y_train2)
test_data2 = GeneDataset2(X_test2, y_test2)

input_size = len(genes_filter_in_adata)
hidden_size1 = 256
hidden_size2 = 128
output_size = 2

# 创建DataLoader
train_loader1 = DataLoader(train_data1, batch_size=32, shuffle=True)
test_loader1 = DataLoader(test_data1, batch_size=32, shuffle=False)

train_loader2 = DataLoader(train_data2, batch_size=32, shuffle=True)
test_loader2 = DataLoader(test_data2, batch_size=32, shuffle=False)

os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575/")

import pickle

# 保存train_loader
with open('train_loader.pkl', 'wb') as f:
    pickle.dump(train_loader, f)

# 保存test_loader
with open('test_loader.pkl', 'wb') as f:
    pickle.dump(test_loader, f)

# 加载train_loader
with open('train_loader.pkl', 'rb') as f:
    train_loader = pickle.load(f)

# 加载test_loader
with open('test_loader.pkl', 'rb') as f:
    test_loader = pickle.load(f)

# 初始化模型和优化器
model1 = DynamicMLP(input_size, hidden_size1, hidden_size2, output_size)
optimizer = torch.optim.Adam(model1.parameters(), lr=0.001)

# 训练模型
for epoch in range(100):
    for i, (inputs, labels) in enumerate(train_loader1):
        outputs = model1(inputs)
        loss = nn.CrossEntropyLoss()(outputs, labels)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
#在黑色素瘤验证集上验证结果
model1.eval()
probs1 = []
true_labels1 = []
with torch.no_grad():
    for inputs, labels in test_loader1:
        outputs = model1(inputs)
        prob = torch.softmax(outputs.data, dim=1)[:, 1]  # 假设正类的标签是1
        probs1.extend(prob.numpy())
        true_labels1.extend(labels.numpy())

from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
# 计算ROC曲线
fpr, tpr, _ = roc_curve(true_labels1, probs1)
roc_auc = auc(fpr, tpr)
#将细胞的概率结果赋值给单细胞数据
adata.obs.cell_probs = probs

# 绘制ROC曲线
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='#A25772', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.legend(prop={'size': 12})
plt.plot([0, 1], [0, 1], color='#7C93C3', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate',fontsize=9)
plt.ylabel('True Positive Rate',fontsize=9)
plt.title('Receiver Operating Characteristic',fontsize=12)
plt.legend(loc="lower right")
plt.show()
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/6.1.Immunotherapy/melanoma_data/GSE120575")

# 保存图片
plt.savefig('GSE120575tumorICI_roc_curve.pdf')

# 保存模型
torch.save({'model_state_dict': model1.state_dict(),
    'optimizer_state_dict': optimizer.state_dict()
}, 'model_checkpoint.pth')

# 加载模型
model = DynamicMLP(input_size, hidden_size1, hidden_size2, output_size)  # 重新创建一个模型实例
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)  # 重新创建一个优化器实例

checkpoint = torch.load('model_checkpoint.pth')
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
####################################
#梯度计算基因权重
# 初始化梯度存储数组
grads = torch.zeros((input_size,))

# 遍历训练数据
for inputs, labels in train_loader:
    inputs.requires_grad = True  # 设置requires_grad=True以计算梯度# 执行前向传播
    outputs = model(inputs)# 选择一个类别计算梯度
    output_idx = outputs.argmax(dim=1)
    output_max = outputs[torch.arange(outputs.shape[0]), output_idx]# 执行反向传播，计算梯度
    output_max.backward(torch.ones_like(output_max))# 累积每个批次的梯度
    grads += inputs.grad.abs().sum(dim=0).detach()# 清空梯度，为下一轮计算做准备
    model.zero_grad()

# 计算平均梯度
feature_importance = grads / len(train_loader.dataset)

# 打印出最重要的10个特征的索引
print(feature_importance.argsort(descending=True)[:10])
important_genes_indices = feature_importance.argsort(descending=True)
important_genes = [genes_filter_in_adata[i] for i in important_genes_indices]
data = {'important_genes': important_genes, 'feature_importance': feature_importance[important_genes_indices]}
df = pd.DataFrame(data)
df.to_csv("feature_importance.csv")

#绘制排序柱状图
import matplotlib.pyplot as plt
plt.clf()
# 取出前10个重要性基因和其对应的重要性
top_10_genes = df.head(10)

# 排序
top_10_genes = top_10_genes.sort_values(by='feature_importance', ascending=True)
plt.figure(figsize=(4, 4))
# 绘制柱状图
plt.barh(top_10_genes['important_genes'], top_10_genes['feature_importance'],color="#9CA777")
plt.xlabel('Feature Importance')
plt.ylabel('Important Genes')
plt.title('Top 10 Important Genes and Their Feature Importance')
plt.grid(False)
plt.tight_layout()
# 保存为 PDF
plt.savefig('top_10_genes_importance.pdf')

# 展示图形
plt.show()
plt.clf()

import matplotlib.pyplot as plt
import numpy as np

# 获取特征重要性排序的索引
sorted_indices = np.argsort(feature_importance)
# 存储每个特征数量下的模型性能
performances = []
best_auc = 0.0
best_model = None

from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score
#你的模型是一个二分类问题，所以你可以使用torch.nn.CrossEntropyLoss作为你的损失函数。你可以在初始化模型和优化器的部分添加如下代码来定义你的损失函数
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
plt.clf()
plt.figure()
plt.plot(range(5, len(sorted_indices) + 1, 5), performances,color="#67729D")
plt.xlabel('Number of Features',fontsize=9)
plt.ylabel('Performance: auc',fontsize=9)
plt.show()
plt.grid(False)
plt.savefig('feature_num_performance.pdf')

# 找出性能最好的特征数量
best_num_features = (np.argmax(performances))
print('Best number of features:', best_num_features)


def MLP_model_test(adata,data_index):
    X_train, X_test, y_train, y_test = train_test_split(adata.X, adata.obs['response_num'], test_size=0.2, random_state=42)
    train_data = GeneDataset(X_train, y_train)
    test_data = GeneDataset(X_test, y_test)
    input_size = len(genes_filter_in_adata)
    hidden_size1 = 256
    hidden_size2 = 128
    output_size = 2
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
    plt.legend(prop={'size': 12})
    plt.plot([0, 1], [0, 1], color='#7C93C3', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate',fontsize=9)
    plt.ylabel('True Positive Rate',fontsize=9)
    plt.title('Receiver Operating Characteristic',fontsize=12)
    plt.legend(loc="lower right")
    plt.savefig(data_index + '_roc_curve.pdf')
    torch.save({data_index + 'model_state_dict': model.state_dict(),'optimizer_state_dict': optimizer.state_dict()}, data_index + 'model_checkpoint.pth')
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
    df.to_csv(data_index + "feature_importance.csv")
    plt.clf()
    top_10_genes = df.head(10)
    top_10_genes = top_10_genes.sort_values(by='feature_importance', ascending=True)
    plt.clf()
    plt.figure(figsize=(4, 4))
    plt.barh(top_10_genes['important_genes'], top_10_genes['feature_importance'],color="#9CA777")
    plt.xlabel('Feature Importance')
    plt.ylabel('Important Genes')
    plt.title('Top 10 Important Genes and Their Feature Importance')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(data_index + 'top_10_genes_importance.pdf')
    return train_loader, test_loader, feature_importance