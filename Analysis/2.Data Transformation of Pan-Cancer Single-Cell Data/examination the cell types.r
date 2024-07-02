library(Seurat)
library(Matrix)
library(scPagwas)

cancers=c('BreastCancer','EndometrialCarcinoma','KidneyCancer','ProstateCancer','ColorectalCancer','LiverCancer','LungCancer','HeadandNeck','Melanoma','OvarianCancer','Pancreatic','EsophagealCancer','BladderCancer','GastricCancer','ThyroidCancer')
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/")
cancers=c('BladderCancer','GastricCancer','ThyroidCancer')
for(i in cancers){
    print(i)
    data_df <- arrow::read_feather(paste0(i,'_matrix.feather'))
    my_r_matrix<-data.matrix(data_df)
    obs<-read.csv(paste0(i,"_obs.csv"),header=T)
    if("merge_celltype_annotation" %in% colnames(obs)){
        print("Merge列存")
    }else{
        print("Merge列不存")

    }
    table(obs$merge_celltype_annotation)
    var<-read.csv(paste0(i,"_var_names.csv"))
    #rownames(obs)<-obs$X
    #colnames(my_r_matrix)<-obs$X
    rownames(my_r_matrix)<-var[,1]
    my_r_matrix <- my_r_matrix[, -ncol(my_r_matrix)]
    obs<-obs[,c("louvain_res1","louvain_1_anno")]
    colnames(obs)[2]<-"cell_type"
    rownames(obs)<-colnames(my_r_matrix)
    #创建seruat格式数据
    data<-CreateSeuratObject(
  my_r_matrix,
  project = "my_r_matrix",
  assay = "RNA",
  names.delim = "_",
  meta.data = obs
)
Idents(data)<-data$merge_celltype_annotation
saveRDS(data,file=paste0(i,'_scdata.rds'))
}

for(i in cancers){
    print(i)
    data<-readRDS(paste0(i,'_scdata.rds'))
    obs<-read.csv(paste0(i,"_obs.csv"),header=T)
    data$merge_celltype_annotation<-obs$merge_celltype_annotation
    Idents(data)<-data$merge_celltype_annotation

    if("merge_celltype_annotation" %in% colnames(data@meta.data)){
        print("Merge列存")
    }
    keywords <- c("B/Plasma", "DC", "Mac/Mono", "Mast", "NK","T.cyt","Th","Treg")
        # 判断 annotation 列是否包含关键字
        is_contained <- unique(Idents(data)) %in% keywords
        # 输出逻辑判断结果
        if(all(is_contained)==TRUE){
             print("Merge列全")
              saveRDS(data,file=paste0(i,'_scdata.rds'))
        }else{
            print("Merge列不全")
        }
}


    saveRDS(data,file=paste0(i,'_scdata.rds'))
