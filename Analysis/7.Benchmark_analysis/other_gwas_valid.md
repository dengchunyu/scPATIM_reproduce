## 利用另一套gwas数据对结果进行验证细胞类型结果

计算得分
```bash
source activate mypy

cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
cancers=('BreastCancer' 'ProstateCancer' 'Melanoma' 'LungCancer' 'BladderCancer' 'KidneyCancer' 'ThyroidCancer' 'OvarianCancer' 'Pancreatic' 'GastricCancer' 'ColorectalCancer' 'EndometrialCarcinoma')
gwass=('ieu-b-4810' 'bbj-a-148' 'ieu-b-4969' 'ieu-b-4954' 'ieu-b-4874' 'ukb-b-1316' 'GCST90018929' 'ieu-b-4963' 'bbj-a-140' 'bbj-a-119' 'ieu-b-4965' 'GCST90018838')
cancer='LungCancer'
gwas='ieu-b-4954'
python scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/${cancer}_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${cancer}/${gwas}_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${cancer}/${gwas}_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${cancer}/${gwas}_celltypeP.csv

```

```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/another_gwas")
BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'

cancers=['BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','ColorectalCancer','EndometrialCarcinoma']
gwass=['ieu-b-4810','bbj-a-148','ieu-b-4969','ieu-b-4954','ieu-b-4874','ukb-b-1316','GCST90018929','ieu-b-4963','bbj-a-140','bbj-a-119','ieu-b-4965','GCST90018838']
def p_merge(pvalues):
    zvalues = -np.sqrt(2) * norm.ppf(pvalues / 2)
    ztotal = np.mean(zvalues)
    p_total = norm.cdf(-abs(ztotal))
    return p_total

def merge_celltype_p(single_p, celltype):
    celltype_p = pd.DataFrame({'celltype': celltype, 'pvalue': single_p})
    celltype_p = celltype_p.groupby('celltype')['pvalue'].agg(p_merge).reset_index()
    return celltype_p
for i in range(5,6):
    print(i)
    gwas=gwass[i]
    cancer=cancers[i]
    results_file=os.path.join(BaseDirectory,cancer+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    drs_df=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/" + cancer+ "/" +gwas+ "_scDRS_score.csv")
    std_multiplier = 2
    drs_df=drs_df.iloc[:,[1,3]]
    drs_df.columns='scDRS_'+drs_df.columns
    drs_df.index=adata.obs.index
    drs_df['first_celltype_annotation']=adata.obs['first_celltype_annotation']
    drs_df['merge_celltype_annotation']=adata.obs["merge_celltype_annotation"]
    obs=drs_df
    df=obs['scDRS_norm_score']
    mean = df.mean()
    std = df.std()
    threshold = mean + std_multiplier * std
    df_cleaned = df.clip(lower=mean - std_multiplier * std, upper=threshold)
    df_cleaned = pd.DataFrame(df_cleaned)
    scaler = StandardScaler()
    obs['scDRS_norm_score'] = scaler.fit_transform(df_cleaned)
    ct_value=merge_celltype_p(obs['scDRS_pval'], adata.obs['first_celltype_annotation'].astype(str))
    ct_value.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/'+cancer+'/'+gwas+'_first_celltype_pvalue.csv')
    ct_value=merge_celltype_p(obs['scDRS_pval'], adata.obs["merge_celltype_annotation"].astype(str))
    ct_value.to_csv('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/'+cancer+'/'+gwas+'_Merged_celltype_pvalue.csv')
    obs['Significant_cells'] = obs['scDRS_pval'].apply(lambda x: 1 if x < 0.05 else 0)
    grouped = obs.groupby('merge_celltype_annotation')['Significant_cells']
    grouped_percentage = (grouped.sum() / grouped.count()) * 100
    average_scores = obs.groupby('merge_celltype_annotation')['scDRS_norm_score'].mean()
    ct_value['sig_percentage']=grouped_percentage[ct_value['celltype']].values
    ct_value['celltype'] = ct_value['celltype'].replace('nan', pd.NA)
    ct_value = ct_value.dropna(subset=['celltype'])
    ct_value['average_trs']=average_scores[ct_value['celltype']].values
    output_file = f'{gwas}_merge_average_pscores.csv'
    ct_value.to_csv(output_file, index=False)

```
可视化结果
```R
library(ggplot2)
p_df<-list()
trs_df<-list()
per_df<-list()
gwass<-c('ieu-b-4810','bbj-a-148','ieu-b-4969','ieu-b-4954','ieu-b-4874','ukb-b-1316','GCST90018929','ieu-b-4963','bbj-a-140','bbj-a-119','ieu-b-4965','GCST90018838')
cancers<-c('BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','ColorectalCancer','EndometrialCarcinoma')
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/another_gwas")
cell_list<-lapply(gwass,function(i){
df<-read.csv(paste0(i,"_merge_average_pscores.csv"))
return(df$celltype)
})

a<-table(unlist(cell_list))
cellnames<-names(a)[a>=7]

data_list<-lapply(gwass,function(i){
df<-read.csv(paste0(i,"_merge_average_pscores.csv"))
df<-df[df$celltype %in% cellnames,]
return(df)
})
names(data_list)<-cancers

p_df<-list()
trs_df<-list()

for(i in cancers){
df<-data_list[[i]]
#df<-read.csv(paste0(i,"_average_pscores.csv"))
pvalue<- -log10(df$pvalue)
p2<-rep(0,length(cellnames))
names(p2)<-cellnames
p2[df$celltype]<-pvalue
p_df[[i]]<-p2
t<- df$average_trs
t2<-rep(0,length(cellnames))
names(t2)<-cellnames
t2[df$celltype]<-t
trs_df[[i]]<-t2
}


p_df<-as.data.frame(p_df)
trs_df<-as.data.frame(trs_df)

#c_df<-as.data.frame(c_df)

save(p_df,trs_df,file="other_gwas_trs_p_result.RData")
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tidyverse)


pmt <- as.matrix(p_df)
col_fun<-colorRamp2(c(-1.2,0,2.7), c("#9CA777","#FEE8B0", "#F97B22"))
mean_trs<-t(as.matrix(trs_df))
pmt<-t(pmt[,rownames(mean_trs)])
pdf("celltype_anothergwas_P_heatmap.pdf",width=6)
Heatmap(mean_trs, col = col_fun,
    row_km =3, column_km =3,
    cell_fun = function(j, i, x, y, w, h, fill) {
    if(pmt[i, j] > 3) {
        grid.text("***", x, y)
    } else if(pmt[i, j] > 2) {
        grid.text("**", x, y)
    }else if(pmt[i, j] > 1.3) {
       grid.text("*", x, y)
}}
)
dev.off()
```

## 对可靠的gwas数据进行相关性绘图

```python
## 不同gwas对结果的验证
cancers=['BreastCancer','LungCancer','KidneyCancer','ColorectalCancer']

gwas1=['finngen_r7_c3_breast_exallc','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_colorectal']
gwas2=['ieu-b-4810','ieu-b-4954','ukb-b-1316','ieu-b-4965']

```

```python
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm
import os
import scipy.sparse
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
from scipy.stats import linregress
from sklearn.preprocessing import StandardScaler

os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/another_gwas")
sc.set_figure_params(fontsize=16, dpi=80, dpi_save=300, format='svg')
for i in range(0,4):
    print(i)
    g1=gwas1[i]
    g2=gwas2[i]
    cancer=cancers[i]
    drs_df1=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/" + cancer+ "/" +g1+ "_scDRS_score.csv")
    drs_df2=pd.read_csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/" + cancer+ "/" +g2+ "_scDRS_score.csv")
    drs_df1=std_score(drs_df1)
    drs_df2=std_score(drs_df2)
    data = pd.DataFrame({g1: drs_df1['norm_score'],g2: drs_df2['norm_score']})
    plt.figure(figsize=(6, 5))
    sns.set_style('whitegrid')
    sns.regplot(x=g1, y=g2, data=data, scatter_kws={'s': 0.5, 'color': '#FCE09B'}, line_kws={'color': '#B2533E'})
    result = linregress(data[g1], data[g2])
    correlation_coefficient = result.rvalue
    p_value = result.pvalue
    equation = f'y = {result.slope:.2f}x + {result.intercept:.2f}'
    plt.text(0.1, 0.9, f'Correlation Coefficient: {correlation_coefficient:.2f}', transform=plt.gca().transAxes)
    plt.text(0.1, 0.8, f'P-Value: {p_value:.4f}', transform=plt.gca().transAxes)
    plt.title(cancer,fontsize=18)
    sns.despine()
    plt.savefig(cancer+'_tow_gwas_cor_plot', dpi=300, bbox_inches='tight')
    plt.close()

def std_score(obs):
    std_multiplier = 4
    df=obs['norm_score']
    mean = df.mean()
    std = df.std()
    threshold = mean + std_multiplier * std
    df_cleaned = df.clip(lower=mean - std_multiplier * std, upper=threshold)
    df_cleaned = pd.DataFrame(df_cleaned)
    scaler = StandardScaler()
    obs['norm_score'] = scaler.fit_transform(df_cleaned)
    return obs
