### 寻找其他例子数据

- 寻找其他的能够作为例子的数据
- 白血病是造血系统中未成熟的免疫细胞（主要是白细胞）在骨髓中异常增殖造成的恶性疾病，常见于淋巴样和髓样系列。比如：
  - **急性淋巴细胞白血病（ALL）**：通常由 B 或 T 淋巴细胞前体细胞异常增殖引起。
  - **急性髓系白血病（AML）**：髓系细胞前体异常增殖。
  - **慢性淋巴细胞白血病（CLL）**：主要涉及 B 细胞异常克隆性增殖。
- 这些细胞原本是负责免疫或造血功能，但由于基因突变失控，导致了癌症。
  
finngen_r7_cd2_follicular_lymphoma_exallc_prune.txt
滤泡性淋巴瘤（follicular lymphoma）通常是由 B淋巴细胞 发生癌变引起的，具体来说，是 生发**中心B细胞（germinal center B cells）的恶性增殖。**这类淋巴瘤的癌细胞形态类似于正常生发中心中的 中心细胞（centrocytes） 和 中心母细胞（centroblasts），并在淋巴结中形成类似滤泡（follicle）的结构。

finngen_r7_cd2_lymphoid_leukaemia_exallc_prune.txt
淋巴细胞白血病（lymphoid leukaemia）是一类白血病，主要涉及**淋巴细胞**。淋巴细胞白血病通常分为急性和慢性两类，急性淋巴细胞白血病（ALL）和慢性淋巴细胞白血病（CLL）是其中最常见的两种类型。

finngen_r7_cd2_primary_lymphoid_hematopoietic_exallc_prune.txt
原发性淋巴造血系统肿瘤（primary lymphoid hematopoietic）是一类罕见的淋巴瘤，通常由淋巴细胞或造血细胞发生癌变引起。原发性淋巴造血系统肿瘤包括多种亚型，如淋巴母细胞瘤（lymphoblastic lymphoma）、淋巴瘤性肉瘤（lymphoma-like myeloma）等。

将这三种肿瘤在BMMC数据中进行计算

## 构造模拟的gwas数据


```R
# 加载必要包
library(dplyr)
library(tidyr)
library(data.table)
library(scPagwas)

eqtl <- fread("/share/pub/dengcy/Cancer_Gwas/Runtime1.0/2.scPagwas_run/2.immune_snp/onek1k_eqtl_filter5e2.tsv") 

library(scPagwas)
load("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/mp_immune_genelst.RData")
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/")
Single_data<-readRDS("Seu_Hema_data.rds")
 #1.start to run the wrapper functions for example.
 Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/finngen_r7_cd2_follicular_lymphoma_exallc_prune.txt", 
                     Single_data =Single_data,
                     output.prefix="", 
                     output.dirs="follicular_lymphoma",
                     block_annotation = block_annotation,# gene position in chromosome is provided by package. default is hg38, block_annotation_hg37 is hg37.
                     assay="RNA", # the assays for scRNA-seq data to use.
                     Pathway_list=mp_immune_genelst,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,
                     singlecell=T, 
                     celltype=F)
save(Pagwas_data, file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/Pagwas_follicular_lymphoma.RData")

#
Pagwas_data<-scPagwas_main2(Pagwas = NULL,
                     gwas_data ="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/finngen_r7_cd2_lymphoid_leukaemia_exallc_prune.txt", 
                     Single_data =Single_data,
                     output.prefix="", 
                     output.dirs="lymphoid_leukaemia",
                     block_annotation = block_annotation,# gene position in chromosome is provided by package. default is hg38, block_annotation_hg37 is hg37.
                     assay="RNA", # the assays for scRNA-seq data to use.
                     Pathway_list=mp_immune_genelst,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,
                     singlecell=T, 
                     celltype=F
)
save(Pagwas_data, file = "/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/Pagwas_lymphoid_leukaemia.RData")
```
### 读取两个数据的结果并且进行可视化展示

```R
follicular_lymphoma<-read.csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/follicular_lymphoma/_Merged_celltype_pvalue.csv")
lymphoid_leukaemia<-read.csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/lymphoid_leukaemia/_Merged_celltype_pvalue.csv")
#整合两个数据
follicular_lymphoma$cancer_type<-"follicular_lymphoma"
lymphoid_leukaemia$cancer_type<-"lymphoid_leukaemia"
all_data<-rbind(follicular_lymphoma,lymphoid_leukaemia)
all_data$log10_pvalue<- -log10(all_data$pvalue)
##绘制柱状图，横坐标为细胞类型，纵坐标为log10(p值)，颜色根据肿瘤类型不同，分别为follicular_lymphoma和lymphoid_leukaemia
library(ggplot2)
#在y轴加入p小于0.05的标记

pdf("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/model_test/BMMC_lymphoma_test.pdf",width=7,height=5)
ggplot(all_data,aes(x=celltype,y=log10_pvalue,fill=cancer_type))+geom_bar(stat="identity",position="dodge")+theme_minimal()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(title="BMMC in leukaemia",x="celltype",y="-log10(pvalue)")+geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")
dev.off()


```

