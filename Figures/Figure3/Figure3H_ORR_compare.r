
setwd("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/Cancer_Gwas/Runtime_final/Figures/Figure2/")
immunetherapy_orr<-read.csv("immunetherapy_orr2.csv")
library(ggplot2)
library(dplyr)
library("ggsignif")
# 画小提琴图
p <- ggplot(immunetherapy_orr, aes(x = Type, y = ORR, fill = Type)) +
  geom_violin() + 
  #geom_jitter(alpha=0.2,size=1)+
  geom_boxplot(width = 0.1) +
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3) +
  theme_minimal()+
  scale_fill_manual(values=c("#DF826C","#9FBB73","#7C93C3"))+
  geom_signif(
      step_increase = c(0.1), #显著性显示从头由低到高排列
                comparisons=list(c("Module1","Module3"),c("Module1","Module2")),
      map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05), #展示显著性
      textsize=7,tip_length=0.03)
ggsave(p,file="Figure3H_ORR_diff_module.pdf",width=6,height=5)
# 进行两两之间的差异分析
pairwise_comparisons <- pairwise.wilcox.test(immunetherapy_orr$ORR, immunetherapy_orr$Type, p.adjust.method = "BH")

# 提取p值
p_values <- pairwise_comparisons$p.value
#0.038
