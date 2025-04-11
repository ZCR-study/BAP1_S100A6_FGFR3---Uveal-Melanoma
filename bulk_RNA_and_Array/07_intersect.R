###########Intersect#############################
##Environment setup and clearing
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./bulk_RNA_and_Array")

# --- Load Required Libraries ---
library(ggVennDiagram) 
library(writexl)        
library(readxl)         
library(ggplot2)       
library(ggsci)          
library(sf)             

# --- Load and Filter Upregulated Genes from Each Dataset ---

# 1. TCGA Mutation Comparison (BAP1 WT vs MUT)
load("./processed_data/GSE176345/TCGA_UM_BAP1_WT_vs_MUT_limma_DEG.Rdata")
TCGA_Mut <- DEG[DEG$change=="UP",];rm(DEG)

# 2. GSE149920 Dataset
load("./processed_data/GSE149920/GSE149920_deg.Rdata")
GSE149920 <- deg[deg$change=="up",];rm(deg)

# 3. GSE22138 Dataset
load("./processed_data/GSE22138/GSE22138_deg.Rdata")
GSE22138 <- deg[deg$change=="up",];rm(deg)

# 4. GSE176345 Dataset
load("./processed_data/GSE176345/GSE176345_deg.Rdata")
GSE176345 <- DEG[DEG$change=="UP",];rm(DEG)

# 5. TCGA Stage Comparison (Early vs Advance)
load("./processed_data/TCGA/TCGA_UM_Early_vs_Advance_limma_DEG.Rdata")
TCGA_Stage <- DEG[DEG$change=="UP",];rm(DEG)

# --- Prepare List of Upregulated Gene Sets ---
df <- list(
  TCGA_Mut = TCGA_Mut$SYMBOL,
  GSE149920 = GSE149920$SYMBOL,
  GSE22138 = GSE22138$probe_id,
  GSE176345 = GSE176345$SYMBOL,
  TCGA_Stage = TCGA_Stage$SYMBOL
)

# --- Find and Display Common Genes Across All Sets ---
common_genes <- Reduce(intersect,df)
print("Genes upregulated in ALL datasets:")
print(common_genes) 


# --- Generate and Save Venn Diagram ---
pdf("./Rplot/Vene.pdf",width = 8,height = 7)
ggVennDiagram(df[1:5],set_size = 8, label_alpha=0,label_size =6,edge_size = 1.2,label ="count") +
  scale_color_brewer(palette = "Paired")+
  scale_fill_gradient(low="white",high = "white",guide="none")
dev.off()