# 对数据噪声的评估结果

## 构造gwas数据和单细胞数据的噪声

```R
library(Seurat)   # 用于处理单细胞数据
library(dplyr)     # 用于数据处理
library(ggplot2)   # 用于可视化
library(Matrix)    # 用于稀疏矩阵处理
library(readr)     # 用于读取数据

# 设定随机种子
set.seed(42)
### 1. 加载数据
# 加载GWAS数据
gwas_data <- read_tsv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt")

#gwas_ad <- read_tsv("/share/pub/zhouyj/brain_TF/GWAS/PD_pagwas.txt")

# 加载单细胞数据
sc_data <- readRDS("/share/pub/dengcy/GWAS_Multiomics/NM_Healthy_pbmc.rds")
set.seed(123)  # 设置随机种子
sampled_rows <- sample(1:nrow(sc_data), size = 10000, replace = FALSE)

# 提取子集数据
sc_data_subset <- sc_data[sampled_rows, ]
table(Idents(sc_data_subset))
saveRDS(sc_data_subset,file="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/pbmc_sc_data_subset.rds")

sc_data<-readRDS("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/pbmc_sc_data_subset.rds")
# 查看数据结构
head(gwas_data)
print(sc_data_subset)

### 2. 添加噪声的函数
add_noise_to_gwas <- function(gwas_df, noise_level) {
  # 为GWAS数据的beta列添加噪声
  beta_col <- gwas_df$beta
  noise <- rnorm(length(beta_col), mean = 0, sd = noise_level * sd(beta_col))
  gwas_df$beta_noisy <- beta_col + noise
  return(gwas_df)
}

library(Matrix)

add_noise_to_sc_data <- function(sc_data, noise_level) {
  # 获取基因表达矩阵
  sc_expr <- GetAssayData(sc_data, slot = "data")
  
  # 确保矩阵为稀疏矩阵（dgCMatrix）
  if (!inherits(sc_expr, "dgCMatrix")) {
    sc_expr <- as(sc_expr, "dgCMatrix")
  }
  
  # 添加噪声：仅对非零元素添加噪声
  non_zero_indices <- which(sc_expr != 0, arr.ind = TRUE)  # 获取非零元素的索引
  noise <- rnorm(nrow(non_zero_indices), mean = 0, sd = noise_level * mean(sc_expr@x))
  
  # 对非零元素添加噪声，并确保非负值
  sc_expr@x <- sc_expr@x + noise  # 直接对稀疏矩阵的非零值向量操作
  sc_expr@x[sc_expr@x < 0] <- 0   # 将负值置为0
  
  # 将带噪声的矩阵赋值回 Seurat 对象
  sc_data_noisy <- sc_data
  sc_data_noisy@assays$RNA@data <- sc_expr
  
  return(sc_data_noisy)
}



setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/")
### 4. 分别添加噪声并测试影响
gwas_noise_levels <- seq(0.002, 0.01, by = 0.002)  # 更小的噪声比例范围
gwas_noise_levels <- seq(0.02, 0.1, by = 0.02)  # 更小的噪声比例范围

# 1. 测试GWAS数据噪声的影响
for (noise_level in gwas_noise_levels) {
  #为GWAS数据添加噪声
  noisy_gwas <- add_noise_to_gwas(gwas_data, noise_level)
  write.csv(noisy_gwas,file=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/gwas_noise",noise_level,".csv"))
}

sc_noise_levels <- c(seq(0.002, 0.01, by = 0.002),seq(0.02, 0.1, by = 0.02),seq(0.2,0.5, by = 0.1))   # 更小的噪声比例范围

for (noise_level in sc_noise_levels) {
  noisy_sc_data <- add_noise_to_sc_data(sc_data, noise_level)
  saveRDS(noisy_sc_data,file=paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scdata_noise",noise_level,".rds"))
}
```

## 分别运行scPtim

原始结果

```R
library(Seurat)   # 用于处理单细胞数据
library(dplyr) 
library(scPagwas)
# 加载必要的库
library(Seurat)
library(data.table)

noisy_gwas <- read.table('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt',sep='\t',head=TRUE)

load("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/mp_immune_genelst.RData")
Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                     gwas_data =noisy_gwas,
                     Single_data ="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/pbmc_sc_data_subset.rds",
                     output.prefix="", 
                     output.dirs="scPatim_orignal",
                     block_annotation = block_annotation,
                     assay="RNA", 
                     Pathway_list=mp_immune_genelst,
                     n.cores=1,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,
                     singlecell=T, 
                     celltype=F)
```

对gwas和单细胞的噪音数据分别进行计算

```R
library(Seurat)   # 用于处理单细胞数据
library(dplyr) 
library(scPagwas)
library(Seurat)
library(data.table)
Args <- commandArgs(T)
gwas_level = print(Args[1])

gwas_level<-as.numeric(gwas_level)

noisy_gwas <- read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/gwas_noise",gwas_level,".csv"))
noisy_gwas$beta<-noisy_gwas$beta_noisy

load("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/mp_immune_genelst.RData")
Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                     gwas_data =noisy_gwas,
                     Single_data ="/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scdata_noise/pbmc_sc_data_subset.rds",
                     output.prefix="", 
                     output.dirs=paste0("scPatim_gwas",gwas_level),
                     block_annotation = block_annotation,
                     assay="RNA", 
                     Pathway_list=mp_immune_genelst,
                     n.cores=1,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,
                     singlecell=T, 
                     celltype=F)


library(Seurat)   # 用于处理单细胞数据
library(dplyr) 
library(scPagwas)
library(Seurat)
library(data.table)
library(readr)
Args <- commandArgs(T)
scdata_level = print(Args[1])

scdata_level<-as.numeric(scdata_level)

gwas_data <- read_tsv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.1.valid_result/DBCLC_result/DLBCL_eqtls_gwas_data.txt")
gwas_data<-as.data.frame(gwas_data)
load("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/mp_immune_genelst.RData")
Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                     gwas_data =gwas_data,
                     Single_data =paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scdata_noise",scdata_level,".rds"),
                     output.prefix="", 
                     output.dirs=paste0("scPatim_scdata",scdata_level),
                     block_annotation = block_annotation,
                     assay="RNA", 
                     Pathway_list=mp_immune_genelst,
                     n.cores=1,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,
                     singlecell=T, 
                     celltype=F)
```


linux环境运行
```shell
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/
items=(0.002 0.004 0.006 0.008 0.01)
items=(0.02 0.04 0.06 0.08 0.1)

for item in "${items[@]}"
do
echo $item
    filename="gwas_noise_${item}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e gwas_$item.err
#SBATCH -o gwas_$item.out
#SBATCH -J gwas_$item
#SBATCH -w in005
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/
Rscript scPatim_gwas_noise.r ${item}
EOF
    sbatch $filename
done


############单细胞数据的循环
items=(0.02 0.04 0.06 0.08 0.1)
items=(0.002 0.004 0.006 0.008 0.01)
items=(0.2 0.3 0.4 0.5)


for item in "${items[@]}"
do
echo $item
    filename="scdata_noise_${item}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e scdata_$item.err
#SBATCH -o scdata_$item.out
#SBATCH -J scdata_$item
#SBATCH -w in006
#SBATCH --mem=100000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate myr
cd /share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/
Rscript scPatim_scdata_noise.r ${item}
EOF
    sbatch $filename
done

```

## 整合并且评估结果

```R

evaluate_results <- function(original_genes, noisy_genes) {
  # 计算基因重叠率 (Jaccard Index)
  intersection <- length(intersect(original_genes, noisy_genes))
  union <- length(union(original_genes, noisy_genes))
  
  jaccard_index <- intersection / union
  return(jaccard_index)
}
original_result<-read.csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scPatim_orignal/_gene_PCC.csv")
original_genes <- original_result$X[order(original_result$weight_pcc, decreasing = T)[1:500]]

### 4. 分别添加噪声并测试影响
#noise_levels <- c(seq(0.002, 0.01, by = 0.002),seq(0.02, 0.1, by = 0.02))  # 更小的噪声比例范围

gwas_noise_levels <- c(seq(0.002, 0.01, by = 0.002),seq(0.02, 0.1, by = 0.02),seq(0.2,0.5, by = 0.1))   # 更小的噪声比例范围
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/")
# 1. 测试GWAS数据噪声的影响
jaccard_results_gwas <- c()
for (noise_level in gwas_noise_levels) {
  # 运行带噪声数据的scPATIM分析
  noisy_result <- read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scPatim_gwas",noise_level,"/_gene_PCC.csv"))
  noisy_genes<-noisy_result$X[order(noisy_result$weight_pcc, decreasing = T)[1:500]]
  # 计算与原始结果的一致性
  jaccard_index <- evaluate_results(original_genes, noisy_genes)
  jaccard_results_gwas <- c(jaccard_results_gwas, jaccard_index)
}

# 2. 测试单细胞数据噪声的影响
sc_noise_levels <- c(seq(0.002, 0.01, by = 0.002),seq(0.02, 0.1, by = 0.02),seq(0.2,0.5, by = 0.1))   # 更小的噪声比例范围

jaccard_results_sc <- c()
for (noise_level in sc_noise_levels) {
    noisy_result <- read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/review_test/noisy_test/scPatim_scdata",noise_level,"/_gene_PCC.csv"))
  noisy_genes<-noisy_result$X[order(noisy_result$weight_pcc, decreasing = T)[1:500]]
  # 计算与原始结果的一致性
  jaccard_index <- evaluate_results(original_genes, noisy_genes)
  jaccard_results_sc <- c(jaccard_results_sc, jaccard_index)
}

### 5. 可视化结果
# 绘制GWAS数据噪声的影响
jaccard_df_gwas <- data.frame(
  noise_level = gwas_noise_levels,
  jaccard_index = jaccard_results_gwas
)
save(jaccard_df_gwas,file="jaccard_df_gwas.RData")
library(ggplot2)
pdf("gwas_noise_jaccard.pdf",width=6,height=6)
ggplot(jaccard_df_gwas, aes(x = noise_level, y = jaccard_index)) +
  geom_line() +
  geom_point() +
  labs(x = "GWAS Data Noise Level", y = "Jaccard Index",
       title = "Impact of GWAS Data Noise on scPATIM Results") +
  theme_minimal()
dev.off()

# 绘制单细胞数据噪声的影响
jaccard_df_sc <- data.frame(
  noise_level = sc_noise_levels,
  jaccard_index = jaccard_results_sc
)
save(jaccard_df_sc,file="jaccard_df_scdata.RData")

pdf("scdata_noise_jaccard.pdf",width=6,height=6)
ggplot(jaccard_df_sc, aes(x = noise_level, y = jaccard_index)) +
  geom_line() +
  geom_point() +
  labs(x = "Single-Cell Data Noise Level", y = "Jaccard Index",
       title = "Impact of Single-Cell Data Noise on scPATIM Results") +
  theme_minimal()
dev.off()

```