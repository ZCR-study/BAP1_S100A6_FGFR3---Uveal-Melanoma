###########GSE149920_data_process#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA_and_Array")

# --- Load Required Libraries ---
library(stringr)
library(limma)
library(edgeR)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Load Processed Data ---
load("./processed_data/GSE149920/GSE1449920_exp_pd.Rdata")

# --- Prepare Data for Limma ---
pd <- pd[c(1:3,13:15),]
exp <- exp[,pd$title]
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 2), ]

group_list=ifelse(str_detect(pd$title,"MM66"),"BAP1_WT","BAP1_MUT")
group_list = factor(group_list,
                    levels = c("BAP1_WT","BAP1_MUT"))


# --- Limma-Voom Analysis Pipeline ---
design <- model.matrix(~0+group_list) 
colnames(design)=levels(group_list)
rownames(design)=colnames(exp)

dge <- DGEList(counts=exp)
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group_list)),collapse = "-") 
constrasts

cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
deg <- DEG


deg <- mutate(deg,probe_id=rownames(deg))
deg$probe_id <- str_sub(deg$probe_id,1,15)
head(deg)
name <- bitr(deg$probe_id,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
deg <- inner_join(deg,name,by=c("probe_id"="ENSEMBL"))
head(deg)


deg <- deg[!duplicated(deg$SYMBOL),]


logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
table(deg$change)

save(deg,file = "./processed_data/GSE149920/GSE149920_deg.Rdata")