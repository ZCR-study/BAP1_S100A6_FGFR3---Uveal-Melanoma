###########Gene Expression Based Survival Analysis (TCGA Data)#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA_and_Array")

# --- Load Data ---
load("./processed_data/TCGA/TCGA_UM_survival_TPM.Rdata")

# --- Load Required Libraries ---
library(survival)
library(survminer)

# --- Data Preparation ---
exprSet <- log2(exprSet+1)
g = "BAP1" #S100A6

meta$gene = ifelse(as.numeric(exprSet[g,]) > median(as.numeric(exprSet[g,])),'high','low')
meta$gene <- factor(meta$gene,levels = c("low","high"))
table(meta$gene)

# --- Kaplan-Meier (KM) Survival Analysis ---
sfit_K = survfit(Surv(time, event)~gene, data=meta)

# --- Generate Kaplan-Meier Plot ---
ggsurvplot(sfit_K,
           data = meta, 
           legend.title = "Gene Expression",
           pval =TRUE, 
           conf.int = TRUE,
           risk.table = F, # Add risk table
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = c("#eb5354", "#1c8041"),
           ggtheme = theme_bw()
)

ggsave(filename = "./Rplot/BAP1.pdf" #"./Rplot/S100A6.pdf"
       ,width = 3.5,height = 3.5) 

dev.off()

# --- Generate Cumulative Hazard Plot ---
ggsurvplot(sfit_K,
           fun = "event", 
           title = NULL,
           risk.table = F, 
           # risk.table.col = "strata"
           pval = TRUE, # p值
           pval.method=TRUE,
           pval.coord = c(52,1), 
           pval.method.coord= c(25,1),
           censor.shape = "", 
           surv.scale = "percent", 
           #legend =c(0.9, 0.1), 
           legend.title ="Group",
           legend="top",
           legend.labs = c("low","high"),
           xlim = c(0, 85), 
           break.time.by = 20, 
           ylim = c(0, 1), 
           palette = c( "#1c8041","#eb5354"), 
           data = meta,
           conf.int = T,
           ggtheme = theme_bw()
)

ggsave(filename = "./Rplot/BAP1_cumulative.pdf",width = 3.5,height = 3.5)