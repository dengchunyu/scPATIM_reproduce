## 获得所有癌症去批次之前的数据
为了获得去批次之前的数据并且保留注释信息，需要重新合并数据获得adata_concat并且倒入去批次之后的数据adata并且做细胞的交集，将adata中的first_celltype_annotation注释的列合并到adata_concat中，写出实现的python代码
items=('EsophagealCancer' 'LungCancer' 'Pancreatic' 'BreastCancer' 'ThyroidCancer' 'BladderCancer' 'ColorectalCancer' 'LiverCancer' 'OvarianCancer' 'EndometrialCarcinoma' 'KidneyCancer' 'GastricCancer'  'ProstateCancer' 'HeadandNeck' 'Melanoma')

```python
import scanpy as sc
import pandas as pd
import os
import numpy as np
import warnings
import gc
import anndata as ad

#breast
cancertype='BreastCancer'
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE176078/adata.h5ad")
adata1.obs['samples']= adata1.obs['orig.ident']
adata1.obs['patients']= adata1.obs['orig.ident']
adata1.obs=adata1.obs.drop(labels=['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mito', 'subtype', 'celltype_subset', 'celltype_minor',],axis=1)
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/GSE161529/adata.h5ad")
adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BreastCancer/E-MTAB-8107/adata.h5ad")
adata4.obs_names_make_unique()
adata2.obs_names_make_unique()

# 1. 重新合并数据，获得去批次之前的数据
adata_concat = adata1.concatenate(adata2, adata4, batch_categories=['GSE176078', 'GSE161529', 'E-MTAB-8107'])
nobatch_data(adata_concat,cancertype='BreastCancer')

#结直肠癌
adata1 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/E-MTAB-8107/adata.h5ad')
adata2 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE164522/adata.h5ad')
adata2 = adata2[adata2.obs.celltype_major != 'CD45-']
adata3 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE144735/adata.h5ad')
adata4 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132257/adata.h5ad')
adata5 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ColorectalCancer/GSE132465/adata.h5ad')
adata_concat = adata1.concatenate(adata2,adata3,adata4,adata5, batch_categories=['GSE164522','E-MTAB-8107', 'GSE144735','GSE132257','GSE132465']) #合并数据集
cancertype='ColorectalCancer'
nobatch_data(adata_concat,cancertype)

#EsophagealCancer
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE156728/adata.h5ad")
adata1.obs['patients']=adata1.obs['patient']
adata1.obs['samples']=adata1.obs['patient']
adata1.obs=adata1.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE145370/adata.h5ad")
adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EsophagealCancer/GSE160269/adata.h5ad")
adata_concat = adata1.concatenate(adata2,adata4, batch_categories=['GSE156728','GSE145370','GSE160269']) #合并数据集
cancertype='EsophagealCancer'
nobatch_data(adata_concat,cancertype)

#LungCancer
adata1 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/GSE131907/adata.h5ad')
adata2 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/Leader_et_al/adata_all.h5ad')
#adata2是我们需要计算的
adata3 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653/adata.h5ad')
adata2.var_names=adata2.var.gene_ids
adata2.obs.patients
ps=[]
for x in adata2.obs.patients:
    x="Patient_"+str(x)
    ps.append(x)
adata2.obs.patients= ps
adata2.obs.samples= ps
adata_concat = adata1.concatenate(adata2,adata3, batch_categories=['GSE131907', 'Leader_et_al','E-MTAB-6653'])
cancertype='LungCancer'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)


#Pancreatic
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/SCP1644/adata.h5ad")
sc.tl.pca(adata1)
adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE155698/adata.h5ad")
ps=[]
for x in adata3.obs.patients:
    x="PancreasPatient_"+str(x)
    ps.append(x)

adata3.obs.patients= ps
adata3.obs.samples= ps

adata4 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE154778/adata.h5ad")

adata5 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Pancreatic/Pancreatic_GSE156728/adata.h5ad")

adata5.obs['patients']=adata5.obs['patient']
adata5.obs['samples']=adata5.obs['patient']
adata5.obs=adata5.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)

adata3.obs_names_make_unique()
adata4.obs_names_make_unique()
adata5.obs_names_make_unique()

adata_concat = adata1.concatenate(adata3,adata4,adata5, batch_categories=['SCP1644', 'GSE155698','GSE154778','GSE156728']) #合并数据集
cancertype='Pancreatic'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)

#BladderCancer
adata1= sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BladderCancer/GSE211388/adata.h5ad")
samples = [name.split('_')[0]+ "_"+ name.split('_')[1] for name in adata1.obs_names]
samples_df = pd.DataFrame({'samples': samples})
value_counts = samples_df['samples'].value_counts()
adata2= sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/BladderCancer/GSE135337/adata.h5ad")
adata_concat = adata1.concatenate(adata2, batch_categories=['GSE211388', 'GSE135337'])
cancertype='BladderCancer'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)

#LiverCancer
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE151530/adata.h5ad")
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LiverCancer/GSE140228/adata.h5ad")
adata1.obs_names_make_unique()
adata2.obs_names_make_unique()
adata1.obs_names = adata1.obs_names + "_GSE151530"
# 将 var_names 转为普通字符串类型
adata2.var_names = adata2.var_names.astype(str)
# 然后再调用 var_names_make_unique()
adata2.var_names_make_unique()
# 使用 anndata.concat 进行合并
adata_concat = adata1.concatenate(adata2,batch_categories=["GSE151530","GSE140228"],index_unique='raise') #合并数据集
cancertype='LiverCancer'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)



#OvarianCancer
adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/E-MTAB-8107/adata.h5ad')
adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/OvarianCancer/OV_FTC_GSE154763/adata.h5ad')

adata2=adata2[adata2.obs['tissue']=='T']
adata2.obs['patients']=adata2.obs['patient']
adata2.obs['samples']=adata2.obs['patient']
adata2.obs=adata2.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch','library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster','source', 'tech', 'UMAP1', 'UMAP2'],axis=1)
adata1.obs=adata1.obs.drop(labels=['batch'],axis=1)
adata_concat = adata1.concatenate(adata2,batch_categories=['E-MTAB-8107', 'GSE154763']) #合并数据集
cancertype='OvarianCancer'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)

#'EndometrialCarcinoma' 
adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma/UCEC_GSE156728/adata.h5ad')
adata1.obs['patients']=adata1.obs['patient']
adata1.obs['samples']=adata1.obs['patient']
adata1.obs=adata1.obs.drop(labels=['cancerType', 'patient', 'libraryID', 'loc', 'meta.cluster', 'platform'],axis=1)
adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/EndometrialCarcinoma/UCEC_GSE154763/adata.h5ad')
adata2.obs['patients']=adata2.obs['patient']
adata2.obs['samples']=adata2.obs['patient']
adata2.obs=adata2.obs.drop(labels=['percent_mito', 'n_counts', 'percent_hsp', 'barcode', 'batch', 'library_id', 'cancer', 'patient', 'tissue', 'n_genes', 'MajorCluster', 'source', 'tech', 'UMAP1', 'UMAP2'],axis=1)
adata_concat = adata1.concatenate(adata2, batch_categories=['GSE156728','GSE154763'])
cancertype='EndometrialCarcinoma'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)


#'KidneyCancer' 
adata1=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/RC_GSE156728/adata.h5ad')
adata2=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/KIDNEY_GSE154763/adata.h5ad')
adata3=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE121638/adata.h5ad')
adata4=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE222703/adata.h5ad')
adata5=sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/KidneyCancer/GSE202813/adata.h5ad')
adata_concat= adata1.concatenate(adata2,adata3,adata4,adata5,batch_categories=["GSE156728","GSE154763","GSE121638","GSE222703","GSE202813"]) #合并数据集
cancertype='KidneyCancer'
adata_concat.obs_names_make_unique()
nobatch_data(adata_concat,cancertype)


#'GastricCancer'
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE206785/adata.h5ad")
for i in ["Fibroblast","Epithelial", "Endothelial","Mural","Glial"]:
	adata1 = adata1[adata1.obs["Type"] != i]
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/GastricCancer/GSE163558/adata.h5ad")
adata_concat = adata1.concatenate(adata2,batch_categories=['GSE206785','GSE163558']) #合并数据集

adata_concat.obs_names_make_unique()
cancertype='GastricCancer'
nobatch_data(adata_concat,cancertype)


#'ProstateCancer' 
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE141445/adata.h5ad")
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE137829/adata.h5ad")
adata3 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE176031/adata.h5ad")
adata4 = sc.read_h5ad('/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/ProstateCancer/GSE237603/adata.h5ad')
adata_concat = adata1.concatenate(adata2,adata3,adata4,batch_categories=['GSE141445', 'GSE137829','GSE176031','GSE237603'])
adata_concat.obs_names_make_unique()
cancertype='ProstateCancer'
nobatch_data(adata_concat,cancertype)


#'HeadandNeck' 
adata1 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE139324/adata.h5ad")
adata2 = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/HeadandNeck/GSE164690/adata.h5ad")
adata_concat = adata1.concatenate(adata2,batch_categories=['GSE139324', 'GSE164690'])
adata_concat.obs_names_make_unique()
cancertype='HeadandNeck'
nobatch_data(adata_concat,cancertype)

#'Melanoma'
adata1=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE123139/adata.h5ad")
adata1.var_names_make_unique()
adata2=sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/Melanoma/GSE115978/adata.h5ad")
adata2.var_names_make_unique()
adata = adata1.concatenate(adata2,batch_categories=["GSE151530","GSE140228"],index_unique='raise') #合并数据集

cancertype='Melanoma'
nobatch_data(adata_concat,cancertype)


def nobatch_data(adata_concat,cancertype):
    results_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/" + cancertype + "_Immune_cell.h5ad"
    adata = sc.read_h5ad(results_file)
    common_cells = adata.obs_names.intersection(adata_concat.obs_names)
    adata = adata[common_cells, :].copy()
    adata_concat = adata_concat[common_cells, :].copy()
    adata_concat.obs['first_celltype_annotation'] = adata.obs['first_celltype_annotation']
    output_file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/" + cancertype + "_nobatch_with_annotation.h5ad"
    adata_concat.write(output_file)

```
## 分析数据整合之前的批次效应影响

```python
import scib
from scib.metrics import silhouette_batch
from sklearn.metrics import silhouette_samples, silhouette_score
import numpy as np
from sklearn.metrics import adjusted_rand_score
import scanpy as sc
import os
import argparse

parser = argparse.ArgumentParser(description='evaluate batch correction')
parser.add_argument('--cancer_type', '-i', type=str, help='cancer type')
args = parser.parse_args()
cancertype = args.cancer_type

sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

def evaluate_batch_correction(adata, batch_key='batch', label_key='first_celltype_annotation', k=50,cancer="breast"):
    """
    Evaluate batch correction effectiveness with visualization and metrics.
    Parameters:
    - adata: AnnData object after batch correction
    - batch_key: str, key in adata.obs indicating batch labels
    - label_key: str, key in adata.obs indicating biological labels (e.g., cell types)
    - k: int, number of neighbors for batch mixing evaluation
    Returns:
    - results: dict, containing evaluation metrics
    """
    results = {}
    # 1. UMAP Visualization
    print("\n1. UMAP Visualization")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, use_rep='X_pca')  # Use PCA representation
    sc.tl.umap(adata)
    # Plot UMAP colored by no batch
    sc.pl.umap(adata, color=batch_key, title="UMAP of Batches After Correction",save=cancer + '_no_batch')
    print("\n3. Silhouette Score for Batch Labels")
    if 'X_pca' in adata.obsm:
        silhouette_avg = silhouette_score(adata.obsm['X_pca'], adata.obs[batch_key])
        results['silhouette_score'] = silhouette_avg
        print(f"Silhouette Score: {silhouette_avg:.3f}")
    else:
        print("PCA representation ('X_pca') not found in AnnData object.")
    # 4. Adjusted Rand Index (ARI)
    print("\n4. Adjusted Rand Index (ARI)")
    if label_key and label_key in adata.obs:
        sc.tl.leiden(adata, resolution=1.0)  # Clustering
        pred_clusters = adata.obs['leiden']
        true_labels = adata.obs[label_key]
        ari_score = adjusted_rand_score(true_labels, pred_clusters)
        results['ari_score'] = ari_score
        print(f"Adjusted Rand Index (ARI): {ari_score:.3f}")
    else:
        print(f"Biological labels ('{label_key}') not provided or not found.")
    print("\nEvaluation Complete.")
    return results

os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/batch_test')

#cancertype='LiverCancer'
adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/" + cancertype + "_nobatch_with_annotation.h5ad")


results = evaluate_batch_correction(adata=adata, batch_key='batch', label_key='first_celltype_annotation',cancer= cancertype, k=50)
with open(cancertype + '_nobatch_output.txt', 'w') as file:
    file.write("item\tscore\n")
    for key, value in results.items():
        file.write(key + '\t' + str(value) + '\n')

```

服务器运行

```shell
items=('EsophagealCancer' 'LungCancer' 'Pancreatic' 'BreastCancer' 'ThyroidCancer' 'BladderCancer' 'ColorectalCancer' 'LiverCancer' 'OvarianCancer' 'EndometrialCarcinoma' 'KidneyCancer' 'GastricCancer'  'ProstateCancer' 'HeadandNeck' 'Melanoma')
nodes=("in007" "in008")
node_count=${#nodes[@]}
counter=0

for item in "${items[@]}"
do
echo $item
    node="${nodes[$counter % $node_count]}"
    filename="batch_${item}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e batch_$item.err
#SBATCH -o batch_$item.out
#SBATCH -J batch_$item
#SBATCH -w $node
#SBATCH --mem=150000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/batch_test
python 2.py --cancer_type ${item}
EOF
    counter=$((counter + 1))
    sbatch $filename
done

for item in "${items[@]}"
do
filename="batch_${item}.sh"
sbatch $filename
done
```

## 分析数据整合之后的批次效应影响

```python
import scib
from scib.metrics import silhouette_batch
from sklearn.metrics import silhouette_samples, silhouette_score
import numpy as np
from sklearn.metrics import adjusted_rand_score
import scanpy as sc
import os
import argparse

parser = argparse.ArgumentParser(description='evaluate batch correction')
parser.add_argument('--cancer_type', '-i', type=str, help='cancer type')
args = parser.parse_args()
cancertype = args.cancer_type

sc.set_figure_params(fontsize=14, dpi=80, dpi_save=300, format='svg')

def evaluate_batch_correction(adata, batch_key='batch', label_key='first_celltype_annotation', k=50,cancer="breast"):
    """
    Evaluate batch correction effectiveness with visualization and metrics.
    Parameters:
    - adata: AnnData object after batch correction
    - batch_key: str, key in adata.obs indicating batch labels
    - label_key: str, key in adata.obs indicating biological labels (e.g., cell types)
    - k: int, number of neighbors for batch mixing evaluation
    Returns:
    - results: dict, containing evaluation metrics
    """
    results = {}
    # 1. UMAP Visualization
    print("\n1. UMAP Visualization")
    #sc.pp.neighbors(adata, use_rep='X_pca')  # Use PCA representation
    #sc.tl.umap(adata)
    # Plot UMAP colored by batch
    sc.pl.umap(adata, color=batch_key, title="UMAP of Batches After Correction",save=cancer + '_batch')
    # 2. Batch Mixing Score
    #print("\n2. Batch Mixing Score")
    #adata.obsm['X_emb'] = adata.obsm['X_pca']  # 使用PCA作为基础数据
    #mixing_score = silhouette_batch(adata, batch_key=batch_key,label_key=batch_key, embed='X_pca')
    #results['batch_mixing_score'] = mixing_score
    #print(f"Batch Mixing Score: {batch_mixing_score:.3f}")
    # 3. Silhouette Score for Batch Labels
    print("\n3. Silhouette Score for Batch Labels")
    if 'X_pca' in adata.obsm:
        silhouette_avg = silhouette_score(adata.obsm['X_pca'], adata.obs[batch_key])
        results['silhouette_score'] = silhouette_avg
        print(f"Silhouette Score: {silhouette_avg:.3f}")
    else:
        print("PCA representation ('X_pca') not found in AnnData object.")
    # 4. Adjusted Rand Index (ARI)
    print("\n4. Adjusted Rand Index (ARI)")
    if label_key and label_key in adata.obs:
        #sc.tl.leiden(adata, resolution=1.0)  # Clustering
        pred_clusters = adata.obs['louvain_res1']
        true_labels = adata.obs[label_key]
        ari_score = adjusted_rand_score(true_labels, pred_clusters)
        results['ari_score'] = ari_score
        print(f"Adjusted Rand Index (ARI): {ari_score:.3f}")
    else:
        print(f"Biological labels ('{label_key}') not provided or not found.")
    print("\nEvaluation Complete.")
    return results

# Example usage:
# results = evaluate_batch_correction(adata_concat, batch_key='batch', label_key='cell_type', k=50)
# print(results)
os.chdir('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/batch_test')

adata = sc.read_h5ad("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/" + cancertype + "_Immune_cell.h5ad")
#sc.external.pp.bbknn(adata, batch_key='datasets')

results = evaluate_batch_correction(adata=adata, batch_key='datasets', label_key='first_celltype_annotation',cancer= cancertype, k=50)
with open(cancertype + 'output.txt', 'w') as file:
    file.write("item\tscore\n")
    for key, value in results.items():
        file.write(key + '\t' + str(value) + '\n')

```

服务器运行

```shell
items=('EsophagealCancer' 'LungCancer' 'Pancreatic' 'BreastCancer' 'ThyroidCancer' 'BladderCancer' 'ColorectalCancer' 'LiverCancer' 'OvarianCancer' 'EndometrialCarcinoma' 'KidneyCancer' 'GastricCancer'  'ProstateCancer')
items=('HeadandNeck' 'Melanoma')
nodes=("in009")
node_count=${#nodes[@]}
counter=0

for item in "${items[@]}"
do
echo $item
    node="${nodes[$counter % $node_count]}"
    filename="batch_${item}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e batch_$item.err
#SBATCH -o batch_$item.out
#SBATCH -J batch_$item
#SBATCH -w $node
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/batch_test
python 1.py --cancer_type ${item}
EOF
    counter=$((counter + 1))
    sbatch $filename
done

```

## 可视化结果

```R
#可视化
# 'ThyroidCancer', 
cancers<-c('EsophagealCancer' ,'LungCancer' ,'Pancreatic', 'BreastCancer','ColorectalCancer', 'LiverCancer' ,'OvarianCancer' ,'KidneyCancer', 'GastricCancer',  'ProstateCancer','HeadandNeck', 'Melanoma')
setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/review_test/batch_test")
result_list<-lapply(cancers,function(i){
    results <- read.table(paste0(i,"output.txt"),header=T)
    results$cancer<-i
    results$type<-"batch"
    return(results)
})
combined_df <- do.call(rbind, result_list)
silhouette_df<-combined_df[combined_df$item=='silhouette_score',]
ari_df<-combined_df[combined_df$item=='ari_score',]


result2_list<-lapply(cancers,function(i){
    results <- read.table(paste0(i,"_nobatch_output.txt"),header=T)
    results$cancer<-i
    results$type<-"no_batch"
    return(results)
})
combined_df <- do.call(rbind, result2_list)
silhouette2_df<-combined_df[combined_df$item=='silhouette_score',]
ari2_df<-combined_df[combined_df$item=='ari_score',]

#silhouette_df<-rbind(silhouette_df,silhouette2_df)

library(ggplot2)
# 绘制折线图
pdf("silhouette_score_batch.pdf",width=7,height=5)
ggplot(silhouette_df, aes(x = cancer, y = score)) +
  geom_line(group = 1, color = "blue") +   # 折线
  geom_point(color = "blue") +            # 数据点
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +  # +0.5 虚线
  labs(
    title = "Silhouette score for batch effect",
    x = "Cancer Type",
    y = "Silhouette Score"
  ) +
  theme_minimal() +  # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转横轴标签
  dev.off()


pdf("ari_score_batch.pdf",width=7,height=5)
ggplot(ari_df, aes(x = cancer, y = score)) +
  geom_line(group = 1, color = "blue") +   # 折线
  geom_point(color = "blue") +            # 数据点
  labs(
    title = "ari score for batch effect",
    x = "Cancer Type",
    y = "ari Score"
  ) +
  theme_minimal() +  # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转横轴标签
  dev.off()

```