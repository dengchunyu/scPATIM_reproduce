#基于显著MP的结果计算每个癌症模块的MP热图并且展示出来
library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(ggthemes)
library(circlize)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
mp_df<- read.csv("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/Cancers_top_HCPs.csv")
module1<-c('LungCancer','EndometrialCarcinoma','GastricCancer','BreastCancer','Melanoma','ThyroidCancer')
module2<-c('LiverCancer','ProstateCancer','HeadandNeck','KidneyCancer','OvarianCancer')
module3<-c('Pancreatic','EsophagealCancer','ColorectalCancer')
n1<-match(module1,cancers)
n2<-match(module2,cancers)
n3<-match(module3,cancers)

cancel_list<-list(n1,n2,n3)
result_list<-list()
for(m in 1:3){
    model<-cancel_list[m]
    mp_list<-list()
    for(i in model){
        cancer<-cancers[i]
        gwas<-gwass[i]
        print(cancer)
        load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/",cancer,"_",gwas,"_sigcell_marker_PaHeri.RData"))
        top5_genes <- sig_marker2 %>%
    top_n(5, wt =avg_log2FC) %>%
    ungroup()
        mp_list[[i]]<-top5_genes$gene
    }

    mp_list<-unique(unlist(mp_list))
    #获得所有癌症类型中的mp，并且整合
    mp_sig_fc <- list()
    max_min_f<-function(Nel_df){
    maxs <- max(Nel_df)
    mins <- min(Nel_df)
    scaled <- as.vector(scale(Nel_df, center = mins, scale = maxs - mins))
    return(scaled)
    }
    for(i in model){
        cancer<-cancers[i]
        gwas<-gwass[i]
        print(cancer)
        load(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/",cancer,"_",gwas,"_sigcell_marker_PaHeri.RData"))
        rownames(sig_marker2)<-sig_marker2$gene
        sig_marker<-sig_marker2[,c("gene","cancers","adjp","avg_log2FC")]
        sig_marker<-sig_marker[mp_list,]
        sig_marker$avg_log2FC<-max_min_f(sig_marker$avg_log2FC)
        mp_sig_fc[[i]]<-sig_marker$avg_log2FC
    }
    #names(mp_sig_fc)<-cancers
    mp_sig_fc <- Filter(Negate(is.null), mp_sig_fc)
    mp_fc_df<- as.data.frame(mp_sig_fc)
    rownames(mp_fc_df)<-mp_list
    colnames(mp_fc_df)<-cancers[sort(n3)]
    mp_fc_df <- as.matrix(mp_fc_df)
    result_list[[m]]<-mp_fc_df
}


setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/mp_analysis/")
col_fun<-colorRamp2(c(0,1), c("#F5E8B7", "#D83F31"))
pdf("Figure3Fa_mp1_cancers_P_heatmap.pdf",width=6,height=7)
Heatmap(result_list[[1]], col = col_fun)
dev.off()

col_fun<-colorRamp2(c(0,1), c("#F5FCCD", "#12486B"))
pdf("Figure3Fb_mp2_cancers_P_heatmap.pdf",width=6,height=7)
Heatmap(result_list[[2]], col = col_fun)
dev.off()

col_fun<-colorRamp2(c(0,1), c("#D0E7D2", "#618264"))
pdf("Figure3Fc_mp3_cancers_P_heatmap.pdf",width=6,height=7)
Heatmap(result_list[[3]], col = col_fun)
dev.off()