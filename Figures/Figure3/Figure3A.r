library(scRNAtoolVis)
cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
gwass<-c('ieu-b-4912','finngen_r7_c3_breast_exallc','finngen_r7_c3_prostate_exallc','finngen_r7_c3_melanoma','finngen_r7_c3_lung_nonsmall_exallc','finngen_r7_c3_kidney_notrenalpelvis_exallc','finngen_r7_c3_thyroid_gland_exallc','finngen_r7_c3_ovary_exallc','finngen_r7_c3_pancreas_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal','bbj-a-117','ebi-a-GCST006464')
gene_list<-list()
for(i in 1:15){
    cancer<-cancers[i]
    gwas<-gwass[i]
  pcc<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",cancer,"/",gwas,"_gene_PCC.csv"))
  pcc$weight_pcc<-pcc$adj_logp * pcc$PCC
  pcc<-pcc[,c("X","weight_pcc","pvalue","adj_pvalue")]
  colnames(pcc)<-c("gene","weight_pcc","p_val","p_val_adj")
pcc$cluster <-cancer
gene_list[[i]]<-pcc
}
gene_df<-Reduce(rbind, gene_list)
cancers<-c('LungCancer','EndometrialCarcinoma','GastricCancer','BreastCancer','Melanoma','ThyroidCancer','Pancreatic','EsophagealCancer','ColorectalCancer','LiverCancer','ProstateCancer','HeadandNeck','KidneyCancer','OvarianCancer')
cancer_color<-data.frame(cancers=cancers,colors=c("#fc8191","#f6ee7d","#8432a5","#da4485","#f03e00","#94e6c1","#ffae73","#004684","#37963c","#e95e50","#ffd125","#87c1f3","#ff714c","#5c33bf","#2fb289"))
gene_df$cluster<-factor(gene_df$cluster,levels=cancers)

gene_df<-gene_df[gene_df$weight_pcc>0,]
gene_df<-gene_df[gene_df$p_val_adj<0.01,]
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis")
save(gene_df,file="gene_df.RData")
load("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/gene_df.RData")

pdf("Figure3A_gene_jjVolcano.pdf",width=17,height=8)
jjVolcano(diffData =gene_df,tile.col=cancer_color$colors,col.type = "adjustP",adjustP.cutoff = 0.001,weight_pcc.cutoff =1,aesCol = cancer_color$colors)
dev.off()