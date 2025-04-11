###########rawdata_process#############################
# --- Environment Setup and Initialization ---
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./single_cell_RNA")

# --- Load Required Libraries ---
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(data.table)

# --- Load GSE176029 Data (UM/LM Samples) and  GSE115469 Data (Normal Liver)---
folders=list.files('./rawdata/GSE176029_UM_LM/',pattern='[12345678]$')
folders

scList = lapply(paste0('./rawdata/GSE176029_UM_LM/',folders),function(folder){ 
  CreateSeuratObject(counts = Read10X(folder), 
                     project = folder,
                     min.cells = 3, min.features = 200)
})

ct <- fread("./rawdata/GSE115469_Normal_Liver/GSE115469_Data.csv",data.table=F)
ct[1:4,1:4]
rownames(ct)=ct[,1]
ct=ct[,-1]
sce =CreateSeuratObject(counts=ct,
                        project="GSE115469",
                        min.cells = 3,
                        min.features = 200 )

# --- Merge Datasets ---
scedata <- merge(scList[[1]], 
                 y = c(scList[[2]],sce),
                 add.cell.ids = c("UM_LM1","UM_LM2","Normal")
)

# --- Quality Control (QC) Metrics Calculation ---
scedata[["percent.mt"]] <- PercentageFeatureSet(scedata,pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
scedata_HB_m <- match(HB.genes, rownames(scedata@assays$RNA))
scedata_HB.genes <- rownames(scedata@assays$RNA)[scedata_HB_m] 
scedata_HB.genes <- HB.genes[!is.na(scedata_HB.genes)] 
scedata[["percent.HB"]]<-PercentageFeatureSet(scedata, features=scedata_HB.genes) 

# --- Cell Filtering ---
scedata <- subset(scedata, subset = nFeature_RNA > 500 & nCount_RNA >1000 & percent.mt < 15 & percent.HB < 5)

# --- Data Preprocessing ---
scedata <- NormalizeData(scedata) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

# --- Batch Correction using Harmony ---
system.time({scedata <- RunHarmony(scedata, group.by.vars = "orig.ident")})

# --- Clustering and Dimensionality Reduction ---
scedata <- FindNeighbors(scedata, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1)
scedata <- RunUMAP(scedata, reduction = "harmony", dims = 1:15)
scedata <- RunTSNE(scedata, reduction = "harmony", dims = 1:15)


# --- Marker Gene Identification ---
DefaultAssay(scedata) <- "RNA"
markers  <- FindAllMarkers(scedata,
                           only.pos = TRUE,
                           min.pct = 0.25, 
                           logfc.threshold = 0.75)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]


# --- Marker Gene Visualization for Annotation ---
library(ggplot2)
marker <- c("DCN","COL1A1","FGF7" ,#Fibroblast
            "ALB","TF","CYP3A4","PCK1",#Hepatocytes
            "EPCAM","KRT19","SPP1","FXYD2",#Cholangiocytes
            "CD3D","CD3E","CD8A",#T cells
            "CD19", "MS4A1", "CD79A",#B Cells
            "IGHG1","MZB1", "SDC1",#Plasma cells
            "CD68","CD14","CD163","AIF1",#Monocytes and Macrophages
            "FGFBP2", "CX3CR1","GNLY",#NK cells
            "VWF", "PECAM1",#Endothelial cells
            "MLANA", "MITF","DCT" #UM tumor cells
)

DotPlot(scedata,features = marker)+coord_flip()


# --- Filter Specific Clusters (Doublets)---
scedata_f <- subset(scedata, subset = seurat_clusters != 23 & seurat_clusters != 20)

scedata_f@meta.data$seurat_clusters_f <-  apply(data.frame(as.numeric(scedata_f@meta.data$seurat_clusters)), 1, function(x){
  ifelse(x>20,x-2,x-1)
} )

scedata_f@meta.data$seurat_clusters_f <- factor(scedata_f$seurat_clusters_f, levels = 0:21)
Idents(scedata_f) <- scedata_f$seurat_clusters_f
scedata <- scedata_f
DotPlot(scedata,features = marker)+coord_flip()


# --- Cell Type Annotation ---
Idents(scedata) <- scedata@meta.data$seurat_clusters_f
new.cluster.ids <- c("0"="T cells",
                     "1"="Hepatocytes",
                     "2"="T cells",
                     "3"="NK cells",
                     "4"="UM tumor cells",
                     "5"="Hepatocytes",
                     "6"="Monocytes and Macrophages",
                     "7"="Endothelial cells",
                     "8"="UM tumor cells", 
                     "9"="Monocytes and Macrophages",
                     "10"="T cells",
                     "11"="Fibroblast cells",
                     "12"="Plasma cells",
                     "13"="UM tumor cells",
                     "14"="Monocytes and Macrophages",
                     "15"="B cells",
                     "16"="Fibroblast cells",
                     "17"="Cholangiocytes",
                     "18"="Hepatocytes",
                     "19"="Plasma cells",
                     "20"="Monocytes and Macrophages",
                     "21"="UM tumor cells")
scedata <- RenameIdents(scedata, new.cluster.ids)                        
scedata$celltype <- scedata@active.ident

# --- Visualization ---
DimPlot(scedata,reduction = "tsne",group.by = "celltype",label = F)
DimPlot(scedata,reduction = "tsne",label = T)
DimPlot(scedata,reduction = "tsne",group.by = "seurat_clusters",raster = FALSE)
DimPlot(scedata,reduction = "tsne",group.by = "orig.ident",raster = FALSE)
DimPlot(scedata,reduction = "tsne",group.by = "group",cols = c("#6BAD65","#5085A0"),raster = FALSE)

# --- Refine Metadata ---
table(scedata@meta.data$celltype)
scedata@meta.data$orig.ident <- str_sub(scedata$orig.ident,nchar(scedata$orig.ident)-5,nchar(scedata$orig.ident))
scedata@meta.data$group <- ifelse(str_detect(rownames(scedata@meta.data),"UM") ,"UM","normal")
scedata@meta.data$group_celltype <- paste0(scedata$group,"_",scedata$celltype)

save(scedata,markers,file = "./Rdata/GSE176029-GSE115469_UM.Rdata")
