# BAP1_S100A6_FGFR3---Uveal-Melanoma
This repository contains the R code for the bioinformatics analyses presented in the paper 'S100A6 promotes liver metastasis by activating FGFR3 signaling in BAP1 deficient uveal melanoma'.

These scripts perform analyses including single-cell RNA-seq data processing (QC, integration via Harmony, clustering, annotation using Seurat), bulk RNA-seq/Array differential expression
analysis (limma/voom), survival analysis (survival, survminer), and visualization (ggplot2, ggsurvplot, ggVennDiagram) primarily
on TCGA and GEO datasets. 

To use the code, ensure R and necessary packages (e.g., Seurat, limma, dplyr, ggplot2,
survminer, ggVennDiagram) are installed. Place input data files in the designated directories (./rawdata/, ./processed_data/).
Set the correct working directory within each script (./single_cell_RNA or ./bulk_RNA) before execution. Scripts will generate
output plots (./Rplot/) and processed data files (./Rdata/).
