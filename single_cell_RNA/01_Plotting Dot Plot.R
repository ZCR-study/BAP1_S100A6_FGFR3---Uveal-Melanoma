###########Plotting Enhanced Dot Plot#############################
# --- Environment Setup and Initialization ---
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./single_cell_RNA")

# Load the previously processed Seurat object and markers
load("./Rdata/GSE176029-GSE115469_UM.Rdata")
scedata@meta.data$celltype <- factor(scedata@meta.data$celltype,levels = c("Fibroblast cells","Hepatocytes","Cholangiocytes",
                                                                           "T cells","B cells","Plasma cells","Monocytes and Macrophages",
                                                                           "NK cells","Endothelial cells","UM tumor cells"))

DefaultAssay(scedata) <- "RNA"
Idents(scedata) <- scedata@meta.data$celltype

# --- Load Required Libraries ---
library(Seurat)
library(dplyr)
library(cols4all)
library(ggplot2)

# --- Define Marker Genes ---
marker <- c("DCN","COL1A1","FGF7" ,#Fibroblast cells
            "ALB","TF","CYP3A4","PCK1",#Hepatocytes
            "EPCAM","KRT19","SPP1","FXYD2",#Cholangiocytes
            "CD3D","CD3E","CD8A",#T cells
            "CD19", "MS4A1", "CD79A",#B cells
            "IGHG1","MZB1", "SDC1",#Plasma cells
            "CD68","CD14","CD163","AIF1",#Monocytes and Macrophages
            "FGFBP2", "CX3CR1","GNLY",#NK cells
            "VWF", "PECAM1",#Endothelial cells
            "MLANA", "MITF","DCT" #UM tumor cells
)


# --- Generate and Customize Dot Plot ---
p0 <- DotPlot(scedata, 
              features = marker) +
  scale_color_continuous_c4a_seq('ag_grn_yl',reverse = T) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
p0


ggsave(filename = "./Rplot/GSE176029-GSE115469_UM_Dotplot.pdf",width = 12,height = 6.15)

