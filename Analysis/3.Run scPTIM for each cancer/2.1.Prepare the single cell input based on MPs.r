library(Seurat)
library(Matrix)
library(scPagwas)

cancers=c('BreastCancer','EndometrialCarcinoma','KidneyCancer','ProstateCancer','ColorectalCancer','LiverCancer','LungCancer','HeadandNeck','Melanoma','OvarianCancer','Pancreatic','EsophagealCancer','BladderCancer','GastricCancer','ThyroidCancer')
cancers=c('ProstateCancer','ColorectalCancer','LiverCancer','LungCancer','HeadandNeck')
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/")

for(i in cancers){
    print(i)
    data_df <- arrow::read_feather(paste0(i,'_matrix.feather'))
    my_r_matrix<-data.matrix(data_df)
    obs<-read.csv(paste0(i,"_obs.csv"),header=T)
    var<-read.csv(paste0(i,"_var_names.csv"))
    #rownames(obs)<-obs$X
    #colnames(my_r_matrix)<-obs$X
    rownames(my_r_matrix)<-var[,1]
    my_r_matrix <- my_r_matrix[, -ncol(my_r_matrix)]
    obs<-obs[,c("louvain_res1","louvain_1_anno")]
    colnames(obs)[2]<-"cell_type"
    rownames(obs)<-colnames(my_r_matrix)
    Pagwas <- list()
    load("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1.scPagwas_run_mp/mp_immune_genelst.RData")
    Pagwas <- scPagwas::scCount_data_input(Pagwas=Pagwas,count_data=my_r_matrix,meta_data=obs,Pathway_list=mp_immune_genelst,min_clustercells=10)
    out_file=paste0(i,"_single.pagwas.RData")
    Pagwas$VariableFeatures<-rownames(Pagwas$data_mat)
    Pagwas <- Pathway_pcascore_run(
    Pagwas = Pagwas,
    Pathway_list = mp_immune_genelst,
    min.pathway.size = 5,
    max.pathway.size = 300
    )

save(Pagwas,file=out_file)
}
