###########TCGA_Mut_deg_data_process#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA_and_Array")

# --- Load Required Libraries ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(edgeR)
library(limma)

# --- Load Processed Data ---
load("./processed_data/TCGA/TCGA_UM_counts.Rdata")
load("./processed_data/TCGA/TCGA_UM_meta_for_limma.Rdata")

# --- Prepare Data for Limma ---
counts = exp[apply(exp, 1, function(x) sum(x > 0) > 60), ]

meta$limma <- ifelse(meta$mutation=="BAP1 MUT","BAP1_MUT","BAP1_WT")
group_list=meta$limma

group_list= factor(group_list,
                   levels = c("BAP1_WT","BAP1_MUT"))
table(group_list)

# --- Limma-Voom Analysis Pipeline ---
design <- model.matrix(~0+group_list) 
colnames(design)=levels(group_list)
rownames(design)=colnames(counts)

dge <- DGEList(counts=counts) 
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group_list)),collapse = "-") 
constrasts
# [1] "tumor-normal"
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
View(DEG)

logFC_cutoff <- 0.5
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG$SYMBOL <- rownames(DEG)


save(DEG,file = "./processed_data/TCGA/TCGA_UM_BAP1_WT_vs_MUT_limma_DEG.Rdata")