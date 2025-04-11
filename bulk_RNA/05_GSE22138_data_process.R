###########GSE22138_data_process#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA")

# --- Load Required Libraries ---
library(limma)
library(dplyr)
library(tidyverse)
library(dplyr)
library(survival) # Survival Analysis
library(survminer) 

# --- Load Processed Data ---
load("./processed_data/GSE22138/GSE22138_exp_pd.Rdata")

# --- Prepare Data for Differential Expression (Metastasis vs. Non-metastasis) ---
group_list=ifelse(str_detect(pd$`metastasis:ch1`,"no"),'Non-metastasis',"Metastasis")

group_list = factor(group_list,
                    levels = c('Non-metastasis',"Metastasis"))
table(group_list)

# --- Limma Differential Expression Analysis Pipeline ---
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)


deg <- mutate(deg,probe_id=rownames(deg))
head(deg)


logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
table(deg$change)

save(deg,file = "./processed_data/GSE22138/GSE22138_deg.Rdata")



# --- Survival Analysis based on Chromosome 3 Status ---
rm(list = ls());gc()

# --- Reload Data for Survival Analysis ---
load("./processed_data/GSE22138/GSE22138_exp_pd.Rdata")

# --- Prepare Metadata for Survival Analysis ---
pd <- pd[pd$`chromosome 3 status:ch1`!="NA",]
pd <- pd[pd$`chromosome 3 status:ch1`!="partial monosomy",]
meta <- pd
meta$time <- as.numeric(meta$`months to endpoint:ch1`)
meta$event <- ifelse(meta$`metastasis:ch1`=="yes",1,0)
meta$group <- factor(meta$`chromosome 3 status:ch1`,levels =  c("disomy","monosomy"))
colnames(meta)



# --- Fit Kaplan-Meier Survival Model ---
fit <- survfit(Surv(time, event) ~ group, data = meta)

# --- Generate and Save Kaplan-Meier Plot ---
ggsurvplot(fit,
           data = meta, 
           legend.title = "Group",
           pval =TRUE, 
           conf.int = TRUE,
           risk.table = F, # Add risk table
           tables.height = 0.2,
           tables.theme = theme_cleantable(),
           palette = c( "#1c8041","#eb5354"),
           ggtheme = theme_bw()
)

ggsave(filename = "./Rplot/chromosome-GSE22138-km_plot.pdf",width = 3.5,height = 3.5)

dev.off()

# --- Generate and Save Cumulative Hazard/Event Plot ---
ggsurvplot(fit,
           fun = "event", # 累计风险概率
           title = NULL,
           risk.table = F, # 风险表
           # risk.table.col = "strata", # 根据分层更改风险表颜色
           pval = TRUE, # p值
           pval.method=TRUE,
           pval.coord = c(62,1), #p 值位置
           pval.method.coord= c(30,1),
           censor.shape = "", # 删失值标记
           surv.scale = "percent", # 百分比显示生存率
           #legend =c(0.9, 0.1), 
           legend.title ="Group",
           legend="top",
           legend.labs = c("disomy","monosomy"),
           xlim = c(0, 100), # 横坐标轴范围，相当于局部放大
           break.time.by = 20, # 横坐标刻度
           ylim = c(0, 1), # 这里只为演示用，y轴0~200%毫无意义。
           palette = c( "#1c8041","#eb5354"), # nejm配色
           data = meta,
           conf.int = T,
           ggtheme = theme_bw()
)

ggsave(filename = "./Rplot/chromosome-GSE22138-cumulative_plot.pdf",width = 3.5,height = 3.5)