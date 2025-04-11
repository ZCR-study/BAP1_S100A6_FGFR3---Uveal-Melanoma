###########Plotting Customized t-SNE Visualization#############################
# --- Environment Setup and Initialization ---
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("./single_cell_RNA")

# --- Load Data ---
load("./Rdata/GSE176029-GSE115469_UM.Rdata")

tsne = scedata@reductions$tsne@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = scedata@meta.data$group)  

head(tsne)


# --- Define Custom Color Palette ---
allcolour <- c('#d25774','#3ca0cf','#c376a7','#ad98c3','#6d6fa0','#83ab8e','#ece399','#405993','#cc7f73','#df5734')

# --- Create ggplot ---
p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type)) +
  
  geom_point(size = 1 , alpha =1 )  +
  
  scale_color_manual(values = allcolour)

p

p2 <- p  +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))
p2


p3 <- p2 +         
  theme(
    legend.title = element_blank(),
    legend.key=element_rect(fill='white'),
    legend.text = element_text(size=20), 
    legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5)))  
p3


p4 <- p3 + 
  geom_segment(aes(x = min(tsne$tSNE_1) , y = min(tsne$tSNE_2) ,
                   xend = min(tsne$tSNE_1) +10, yend = min(tsne$tSNE_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(tsne$tSNE_1)  , y = min(tsne$tSNE_2)  ,
                   xend = min(tsne$tSNE_1) , yend = min(tsne$tSNE_2) + 10),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(tsne$tSNE_1) +5, y = min(tsne$tSNE_2) -3, label = "tSNE_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(tsne$tSNE_1) -3, y = min(tsne$tSNE_2) + 5, label = "tSNE_2",
           color="black",size = 3, fontface="bold" ,angle=90) 

p4

cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

library(ggrepel)

p4 +
  
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"))

p4 +
  
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none")

p4


ggsave("./Rplot/GSE176029-GSE115469_UM_Dimplot_normal_group.pdf",width = 8,height = 6.2)
