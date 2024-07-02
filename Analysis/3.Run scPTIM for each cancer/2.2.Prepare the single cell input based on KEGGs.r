library(Seurat)
library(Matrix)
library(scPagwas)
cancers=c('ProstateCancer','ColorectalCancer','LiverCancer','LungCancer','HeadandNeck','Melanoma','OvarianCancer','Pancreatic','EsophagealCancer','BladderCancer','GastricCancer','ThyroidCancer')
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/")
for(i in cancers){
    data_df <- arrow::read_feather(paste0(i,'_matrix.feather'))
    my_r_matrix<-data.matrix(data_df)
    obs<-read.csv(paste0(i,"_obs.csv"),header=T)
    var<-read.csv(paste0(i,"_var_names.csv"))
    rownames(my_r_matrix)<-var[,1]
    my_r_matrix <- my_r_matrix[, -ncol(my_r_matrix)]
    obs<-obs[,c("louvain_res1","louvain_1_anno")]
    colnames(obs)[2]<-"cell_type"
    rownames(obs)<-colnames(my_r_matrix)
    Pagwas <- list()
    Pagwas <- scPagwas::scCount_data_input(Pagwas=Pagwas,count_data=my_r_matrix,meta_data=obs,Pathway_list=Genes_by_pathway_kegg,min_clustercells=10)
    out_file=paste0('/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/KEGG_scdata/',i,"_single.pagwas.RData")
    #save(Pagwas,file=out_file)
Pagwas$VariableFeatures<-rownames(Pagwas$data_mat)
Pagwas <- Pathway_pcascore_run(
      Pagwas = Pagwas,
      Pathway_list = Genes_by_pathway_kegg,
      min.pathway.size = 5,
      max.pathway.size = 300
    )
save(Pagwas,file=out_file)
}