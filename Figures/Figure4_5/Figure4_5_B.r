library(ggrepel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(ggrepel)
library(ggpubr)
setwd("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
cancers<-c('Pancreatic','BreastCancer','HeadandNeck','EsophagealCancer','ThyroidCancer','BladderCancer','ColorectalCancer','LiverCancer','LungCancer','Melanoma','OvarianCancer','EndometrialCarcinoma','KidneyCancer','ProstateCancer','GastricCancer')
cancers<-c('Melanoma','BreastCancer')
gwass<-c('finngen_r7_c3_melanoma','finngen_r7_c3_breast_exallc')
p_list<-list()
#pcc_list<-list()
for(n in 1:2){
i<-cancers[n]
j<-gwass[n]
  PCC<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",i,"/",j,"_gene_PCC.csv"))
  df<-data.frame(pcc=PCC$PCC,pvalue=PCC$adj_pvalue,weight=PCC$weight_pcc,gene=PCC$X)
  df$Direction = as.factor(ifelse(df$pvalue < 0.05 & df$pcc > 0.1,'SignificantPCC','NOSignificant'))
  df$cancer<-i
  df$label <- NA
  df$label[order(PCC$weight_pcc, decreasing =T)[1:10]]<-df$gene[order(PCC$weight_pcc, decreasing =T)[1:10]]
  p<-ggplot(df, aes(x=pcc, y=-log10(pvalue),color=weight,size=weight)) + 
    geom_point(alpha=0.8) + 
    theme_bw(base_size = 12) + 
    xlab("PCC") +
    ylab("LogP") +
    labs(title=i)+
    theme(plot.title = element_text(size=15,hjust = 0.5))
    p<-p+scale_colour_gradient(low= "#F8F1F1", high="#E57C23")
    
  p<-p+geom_hline(yintercept = -log10(0.05), lty = 4) +
    geom_vline(xintercept = 0.1, lty = 4)+
    geom_label_repel(data = df, aes(label = label),
                      colour ="#454545",
                     size = 5,box.padding = unit(0.5, "lines"),
                     point.padding = unit(0.8, "lines"),
                     segment.color = "black",
                     show.legend = FALSE, max.overlaps = 10000)
p_list[[i]]<-p
}

#pdf_file <- "scPagwas_gene_volcanoplot.tiff"
#tiff(file = pdf_file, width =30, height =16,units="in",res=300)
#grid.arrange(grobs = p_list, ncol =5)
#dev.off()

pdf("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/Figure4B.pdf", width =6, height =5)
print(p_list[[1]])
dev.off()

pdf("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/Figure5B.pdf", width =6, height =5)
print(p_list[[2]])
dev.off()