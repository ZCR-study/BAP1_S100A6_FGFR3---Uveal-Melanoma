###########GSE176345_data_process#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA")

# --- Load Required Libraries ---
library(limma)
library(edgeR)

# --- Load Processed Data ---
load("./processed_data/GSE176345/GSE176345_exp_pd.Rdata")

# --- Prepare Data for Limma ---
group_list=c(rep("Normal",3),rep("Tumor",3))
group_list= factor(group_list,
                   levels = c("Normal", "Tumor" ))
table(group_list)

# --- Limma-Voom Analysis Pipeline ---
design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(counts)

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
View(DEG)

logFC_cutoff <- 0.5
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
DEG$SYMBOL <- rownames(DEG)

save(DEG,file = "./processed_data/GSE176345/GSE176345_deg.Rdata")