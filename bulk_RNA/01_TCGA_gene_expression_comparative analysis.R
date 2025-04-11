###########Comparative analysis of gene expression levels between two groups#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA")

##between "Disomy chr3" group and "Monosomy chr3" group

# --- Load Data ---
load("./processed_data/TCGA/TCGA-UVM_SNP.Rdata")
load("./processed_data/TCGA/TCGA_UM_survival_TPM.Rdata")

# --- Load Required Libraries ---
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# --- Process BAP1 Mutation Data ---
BAP1 <- data[data$Hugo_Symbol=="BAP1",]
BAP1 <- BAP1[,c("Tumor_Sample_Barcode","Hugo_Symbol")]
BAP1$Tumor_Sample_Barcode <- str_sub(BAP1$Tumor_Sample_Barcode,1,12)
names(BAP1)[2] <- "mutation"
BAP1$mutation <- paste0(BAP1$mutation ," MUT")

# --- Process Expression Data ---
exp <- exprSet
exp <- exp[!duplicated(rownames(exp)),]
exp <- na.omit(exp)
exp <- log2(exp+1)

# --- Merge Mutation Status into Metadata ---
meta <- left_join(meta,BAP1,by=c("ID"="Tumor_Sample_Barcode"))
meta$mutation <- ifelse(!is.na(meta$mutation),"BAP1 MUT","BAP1 WT")

# --- Prepare Data for BAP1 WT vs MUT Comparison ---
group_list=meta$mutation
#设置参考水平，对照在前，处理在后
group_list= factor(group_list,
                   levels = c("BAP1 WT","BAP1 MUT"))
table(group_list)

identical(meta$ID,str_sub(colnames(exp),1,12))

S100A6=exp["S100A6",]

log2_TPM_S100A6 <- t(S100A6)

group <- group_list
dd <- data.frame(group,log2_TPM_S100A6)
#dd
#str(dd)
dd$group = as.factor(dd$group)

# --- Plot S100A6 Expression: BAP1 WT vs MUT ---
ggplot(dd,
       mapping=aes(x=group,y=S100A6, fill=group))+
  geom_jitter(width = 0.15, shape=21, size=3, color="#252a32") + 
  geom_boxplot(staplewidth = 0.5, width=0.75, 
               outliers = FALSE, alpha=0.75) +
  scale_fill_manual(values = c("#1c8041","#eb5354"),
                    guide = "none")  +   ylab("log2 TPM S100A6")+
  stat_compare_means(method = "t.test") + 
  scale_y_continuous(limits = c(6,13),breaks = seq(6,12,2),expand = c(0,0.2)) +
  theme_classic2()

ggsave("./Rplot/S100A6_in_BAP1_WT_vs_BAP1_MUT.pdf",width = 3.85,height = 3.67,
       units = "in",dpi = 600)

dev.off()

# --- Prepare Data for Disomy vs Monosomy Comparison ---
meta$chr <- ifelse(str_detect(meta$chromosome,"3"),"Monosomy 3","Disomy 3")
group_list=meta$chr
group_list= factor(group_list,
                   levels = c("Disomy 3","Monosomy 3"))
table(group_list)
group <- group_list
dd <- data.frame(group,log2_TPM_S100A6)
#dd
#str(dd)
dd$group = as.factor(dd$group)

# --- Plot S100A6 Expression: Disomy vs Monosomy ---
ggplot(dd,
       mapping=aes(x=group,y=S100A6, fill=group))+
  geom_jitter(width = 0.15, shape=21, size=3, color="#252a32") + 
  geom_boxplot(staplewidth = 0.5, width=0.75, 
               outliers = FALSE, alpha=0.75) +
  scale_fill_manual(values = c("#1c8041","#eb5354"),
                    guide = "none")  +   ylab("log2 TPM S100A6")+
  stat_compare_means(method = "t.test") + 
  scale_y_continuous(limits = c(6,13),breaks = seq(6,12,2),expand = c(0,0.2)) +
  theme_classic2()

ggsave("./Rplot/S100A6_in_Disomy chr3_vs_Monosomy chr3.pdf",width = 3.85,height = 3.67,
       units = "in",dpi = 600)
