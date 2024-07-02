cancers<-c('HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma')
library(ggplot2)
p_df<-list()
trs_df<-list()
per_df<-list()
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
cell_list<-lapply(cancers,function(i){
df<-read.csv(paste0(i,"_merge_average_pscores.csv"))
return(df$celltype)
})

a<-table(unlist(cell_list))
cellnames<-names(a)[a>=10]

data_list<-lapply(cancers,function(i){
df<-read.csv(paste0(i,"_merge_average_pscores.csv"))
df<-df[df$celltype %in% cellnames,]
return(df)
})
names(data_list)<-cancers

p_df<-list()
trs_df<-list()
per_df<-list()

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
s<-df$sig_percentage
s2<-rep(0,length(cellnames))
names(s2)<-cellnames
s2[df$celltype]<-s
per_df[[i]]<-s2
}


p_df<-as.data.frame(p_df)
trs_df<-as.data.frame(trs_df)
per_df<-as.data.frame(per_df)
#c_df<-as.data.frame(c_df)

save(p_df,trs_df,per_df,file="trs_p_result.RData")

load("D:/Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/trs_p_result.RData")
setwd("D:/Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tidyverse)

pmt <- as.matrix(p_df)
col_fun<-colorRamp2(c(-1,0,3), c("#9CA777","#FEE8B0", "#F97B22"))
mean_trs<-t(as.matrix(trs_df))
pmt<-t(pmt[,rownames(mean_trs)])
pdf("celltype_scPagwas_P_heatmap.pdf",width=6)
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