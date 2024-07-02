
setwd("D:\\OneDrive\\Cancer_Gwas\\Runtime2.0\\3.scPagwas_analysis\\")
malignant_mp<-read.csv("MP41_50genes.csv")
library(readxl)

# Replace "file.xlsx" with the name of your Excel file
# Replace NULL with the name(s) of the sheet(s) you want to read, or leave NULL to read all sheets
gene_list<-list()
for(i in c("B_cells","Macrophages","CD8","CD4")){
b_mp <- readxl::read_excel("immune_mp_genes.xlsx",sheet =i)
b_mp<-as.data.frame(b_mp)
colnames(b_mp)<-paste0(i,"_",colnames(b_mp))
gene_list[[i]]<-b_mp
}
gene_df<-Reduce(cbind,gene_list)
mp_list<-cbind(malignant_mp,gene_df)
mp_list<-as.list(mp_list)
immune_df<-read.csv("immunecelltype_genelist.csv")
immune_list<-lapply(unique(immune_df$Cell.type),function(i){
	immune_df$Gene[immune_df$Cell.type==i]
})
names(immune_list)<-unique(immune_df$Cell.type)
mp_immune_genelst<-c(mp_list,immune_list)
length(unique(unlist(mp_immune_genelst)))
save(mp_immune_genelst,file="mp_immune_genelst.RData")