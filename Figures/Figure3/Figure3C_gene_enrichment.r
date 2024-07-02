enrich_result_df<-read.csv("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/enrich_result_df.csv")
setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")
select_pa<-c("interleukin-2 production",
"MHC protein binding",
"immunological synapse",
"T cell receptor complex",
"T cell receptor binding",
"MHC class II protein complex",
"antigen binding",
"phagocytic vesicle",
"mast cell mediated immunity",
"respiratory burst",
"chemokine production",
"macrophage activation",
"ficolin-1-rich granule membrane",
"endoplasmic reticulum to cytosol transport",
"protein exit from endoplasmic reticulum",
"DNA-dependent ATPase activity")

heat_gg<- enrich_result_df[enrich_result_df$description %in% select_pa, c("description","logFDR","cancer")]

library(reshape2)
# 假设数据框为heat_gg
#heat_gg
matrix_heatmap <- dcast(heat_gg, description ~ cancer, value.var = "logFDR", fill = 0)
rownames(matrix_heatmap)<-matrix_heatmap$description
matrix_heatmap<-matrix_heatmap[,-1]
matrix_heatmap<-matrix_heatmap[select_pa,]
#matrix_heatmap<-matrix_heatmap[top_df$description,]
library(ComplexHeatmap)
#ha = HeatmapAnnotation(
#    pathway_type =top_df$method ,
#    col = list(Celltype = c(CC="#61A3BA",pos="#EC8F5E",MF="#A2C579")))
matrix_heatmap[matrix_heatmap>4]<-4
matrix_heatmap<-matrix_heatmap[,c("GastricCancer","Melanoma","BreastCancer","LungCancer","ThyroidCancer","EndometrialCarcinoma","LiverCancer","ProstateCancer","KidneyCancer","OvarianCancer","HeadandNeck","EsophagealCancer","ColorectalCancer","Pancreatic")]

pdf("Figure3C_heatmap_allcancer_pathway.pdf",width=7,height=6)
ComplexHeatmap::Heatmap(t(as.matrix(matrix_heatmap)),name ="-log10(FDR)",col=c("#F5F7F8", "#45474B"),border = NA,column_title=NULL,cluster_rows =FALSE,
        cluster_columns =TRUE,
        show_column_dend=TRUE)
dev.off()
