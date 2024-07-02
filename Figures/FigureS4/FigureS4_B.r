library(ggrepel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(ggrepel)
library(ggpubr)
p2_list<-list()
cancers<-c('ThyroidCancer','LungCancer','Melanoma','BreastCancer','GastricCancer','ColorectalCancer','Pancreatic','EsophagealCancer','LiverCancer','KidneyCancer','ProstateCancer','OvarianCancer','EndometrialCarcinoma','HeadandNeck')
cancer_color<-data.frame(cancers=cancers,colors=c("#fc8191","#f6ee7d","#8432a5","#da4485","#f03e00","#94e6c1","#ffae73","#37963c","#e95e50","#ffd125","#87c1f3","#ff714c","#5c33bf","#2fb289"))

for(i in cancers){

    trs_df<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/",i,"/_singlecell_scPagwas_score_pvalue.Result.csv"))

  df<-data.frame(gPas =trs_df[,'scPagwas.gPAS.score'])
    quantile_95 <- quantile(df$gPas, 0.95)
  p<-ggplot(df, aes(x = gPas)) +
    geom_density(fill = cancer_color[cancer_color$cancers == i,"colors"], alpha = 0.7)+ labs(title=i,x="gPas")+ theme_void()+
    geom_vline(xintercept = quantile_95, linetype = "dashed", color = "red", size = 1)
    #scale_x_continuous(limits = lim)
  
  p2_list[[i]]<-p
}
pdf_file <- "gPas_density.pdf"
pdf(file = pdf_file, width =20, height =10)
grid.arrange(grobs = p2_list, ncol =5)
dev.off()