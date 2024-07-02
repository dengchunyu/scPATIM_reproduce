setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/")
mm_df<-read.csv("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime2.0/3.1.valid_result/cancer_mult_method_presult.csv")
head(mm_df)
library(ggplot2)
library(gridExtra)
library(ggthemes)


    dd<-mm_df[mm_df$Cancer==unique(mm_df$Cancer)[1],]
p <- ggplot(dd, aes(x = Celltype, y = logp, fill = Method)) +
  geom_col(position = "dodge", width = 0.7)+
  labs(x = "Cell Type", y = "logp", fill = "Method")+  
  scale_fill_manual(values=c('#B0C5A4','#D37676','#EBC49F','#6196A6','#F6995C','#AD88C6'))+
theme_few() +
    coord_flip()+
       geom_hline(aes(yintercept=1.3), linetype="dashed", color = "black")+
   pdf("test.pdf", width = 5, height =8)
p
dev.off()    


p_list<-list()
for(i in unique(mm_df$Cancer)){
    dd<-mm_df[mm_df$Cancer==i,]
p_list[[i]] <- ggplot(dd, aes(x = Celltype, y = logp, fill = Method)) +
  geom_col(position = "dodge", width = 0.7)+
  labs(title=i,x = "Cell Type", y = "logp", fill = "Method")+  
  scale_fill_manual(values=c('#B0C5A4','#D37676','#EBC49F','#6196A6','#F6995C','#AD88C6'))+
theme_few() +
    coord_flip()+
       geom_hline(aes(yintercept=1.3), linetype="dashed", color = "black")+
       theme(legend.position="none")
}
pdf("cancer_mult_method_presult_barplots.pdf", width = 10, height =12)
gridExtra::grid.arrange(grobs = p_list, ncol = 5)
dev.off()