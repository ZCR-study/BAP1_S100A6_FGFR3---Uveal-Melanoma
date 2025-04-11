###########Plotting Cell Type Ratios Between Groups#############################
# --- Environment Setup and Initialization ---
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./single_cell_RNA")


# --- Environment Setup and Data Loading ---
load("./Rdata/GSE176029-GSE115469_UM.Rdata")

# --- Data Filtering ---
# Remove UM tumor cells to focus on comparing non-tumor cell types between groups
scedata <- subset(scedata,subset = celltype!="UM tumor cells")

# --- Prepare Data for Ratio Calculation ---
the_metatable=scedata@meta.data[,c("group","celltype")]
colnames(the_metatable)=c("group","celltype")

# --- Calculate Cell Counts per Group and Cell Type ---
tmp2 = the_metatable
tmp2 = tmp2[!duplicated(tmp2),]
rownames(tmp2) = NULL

# Initialize a column to store the counts
tmp2$cell_num <- 0
for(t in rownames(tmp2)){  
  tmp2[t,"cell_num"] = nrow(the_metatable[the_metatable$group == tmp2[t,"group"] 
                                          & the_metatable$celltype == tmp2[t,"celltype"],])
}

table(scedata$celltype)

# --- Load Necessary Libraries ---
library(ggalluvial)
library(ggplot2)
library(gridExtra)
library(plyr)

# --- Calculate Percentage Within Each Cell Type ---
df <- tmp2
df$celltype <- factor(df$celltype,levels=rev(c("T cells","Hepatocytes","Monocytes and Macrophages","NK cells",
                                               "Endothelial cells","Fibroblast cells","Plasma cells","B cells","Cholangiocytes")))
df$group <-factor(df$group,levels=c("normal","UM")) 
df2 = ddply(df,'celltype',transform,percent=cell_num/sum(cell_num))  

c("#71C68B","#E28F93")
# --- Create Plots ---
pdf("test.pdf", width=10, height=6)  


p1<-ggplot(df,aes(x=celltype,y=cell_num,fill=group)) +
  geom_bar(stat='identity',position="stack")+
  scale_fill_manual(values = c("#71C68B","#E28F93"))+ 
  labs(x = "", y="Cell Number") +#设置轴标题
  coord_flip()+   #横向图
  theme_classic() +  
  theme(axis.title.x = element_text(size=20, face='bold'), 
        axis.text.x = element_text(size=16),  
        axis.text.y = element_text(size=16, face='bold'),
        axis.line.x = element_line(size=1.2),
        axis.line.y = element_line(size=1.2),    
        legend.position = c(0.60, 0.12), 
        legend.title = element_text(size=16,face='bold'))+ 
  guides(fill = guide_legend(nrow=1, byrow=TRUE)) 



p2<-ggplot(df2,aes(x=celltype,y=percent,fill=group)) +  
  geom_bar(stat='identity',position="stack")+
  scale_fill_manual(values = c("#71C68B","#E28F93"))+ 
  labs(x = "", y="Sample Ratio") +  
  coord_flip()+  
  theme_classic()+  
  theme(axis.title.x = element_text(size=20, face='bold'),    
        axis.text.x =element_text(size=16),    
        axis.line.x = element_line(size=1.2),    
        axis.line.y = element_line(size=1.2),    
        axis.text.y = element_blank(),    
        legend.position="none") 

# --- Arrange and Save Plots ---
grid.arrange(p1, p2,
             ncol = 7, layout_matrix=rbind(c(1,1,1,1,1,2,2), c(1,1,1,1,1,2,2)
             )) 
ggsave("./Rplot/cellratio.pdf")

dev.off()


