#PCC
#for(cancer in cancers){
    #scmat_path <- file.path(glue::glue("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.scPagwas_run/{cancer}_single.pagwas.RData"))
    #load(scmat_path)
    #score_df<-read.csv(paste0("./", cancer, "/_singlecell_scPagwas_score_pvalue.Result.csv"))
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/")
cancers<-c('ProstateCancer','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','ColorectalCancer','EndometrialCarcinoma')
gwass<-c('bbj-a-148','ieu-b-4954','ieu-b-4874','ukb-b-1316','GCST90018929','ieu-b-4963','bbj-a-140','bbj-a-119','ieu-b-4965','GCST90018838')
for(num in 1:14){
    cancer<-cancers[num]
    gwas<-gwass[num]
    print(cancer)
    load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_pagwas.Rdata"))
    #Pagwas$Pathway_single_results 
    
    scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
    scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
    pcc <- Random_PCC(gPas=scPagwas.gPAS.score,datamat=Pagwas$data_mat,seed=1234,random_times=100,select_num=10000)
    colnames(pcc)<-'PCC'
    top5_index<-order(scPagwas.gPAS.score,decreasing=T)[1:(0.05*(length(scPagwas.gPAS.score)))]
    top5_cellnames<-names(scPagwas.gPAS.score)[top5_index]
    top5_index2<-names(scPagwas.gPAS.score) %in% top5_cellnames
    p_list<-apply(Pagwas$data_mat,1,function(x){
        result <- wilcox.test( x[top5_index2], x[!top5_index2], alternative = "greater")
        adjusted_p_value <- p.adjust(result$p.value, method = "bonferroni")
        return(adjusted_p_value)
        })
    pcc<-as.data.frame(pcc)
    pcc$pvalue <- p_list
    pcc$adj_pvalue <- p.adjust(pcc$pvalue, method = "bonferroni")
    pcc$adj_logp<- -log10(pcc$adj_pvalue)
    pcc$adj_logp[!is.finite(pcc$adj_logp)]<- max(pcc$adj_logp[is.finite(pcc$adj_logp)])+1
    pcc$weight_pcc<- pcc$adj_logp * pcc$PCC
    utils::write.csv(pcc,file = paste0("./", cancer, "/",gwas,"_gene_PCC.csv"), quote = F)
########################
#计算得分
    n_topgenes=500
  assay="RNA"
  scPagwas_topgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = T)[1:n_topgenes]]
  scPagwas_downgenes <- rownames(pcc)[order(pcc$weight_pcc, decreasing = F)[1:n_topgenes]]

  Single_data<-CreateSeuratObject(
    Pagwas$data_mat,
    project = "CreateSeuratObject",
    assay = "RNA",
    names.field = 1,
    names.delim = "_"
  )
  Single_data <- Seurat::AddModuleScore(Single_data, assay = assay, list(scPagwas_topgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.downTRS.Score"))
   a <- data.frame(
    scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1,
    scPagwas.downTRS.Score = Single_data$scPagwas.downTRS.Score2,
    scPagwas.gPAS.score = scPagwas.gPAS.score
  )

    utils::write.csv(a,file = paste0("./", cancer, "/",gwas,"_singlecell_scPagwas_score_pvalue.Result.csv"),quote = F)

}

Random_PCC<- function(gPas,datamat,seed=1234,random_times=100,select_num=10000){
    message("Start to run the corSparse function in random!")
    set.seed(seed)
    sparse_cor_list<-list()
    for (j in 1:random_times) {
      print(paste0("Randome Times: ",j))
    index<-sample(1:ncol(datamat),select_num)
    gPas_select <- data.matrix(gPas[index])
    sparse_cor <- scPagwas::corSparse(
      X = t(datamat[,index]),
      Y = gPas_select
    )
    colnames(sparse_cor) <- "PCC"
    sparse_cor[is.nan(sparse_cor)] <- 0
    sparse_cor_list[[j]]<-unlist(sparse_cor[,1])
    }
    sparse_cor_list<-as.data.frame(sparse_cor_list)
    sparse_cor<- apply(sparse_cor_list,1,mean)
    return(data.matrix(sparse_cor))
}