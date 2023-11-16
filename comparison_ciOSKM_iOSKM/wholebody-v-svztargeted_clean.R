# Compare enriched pathways from GSEA for whole-body reprogramming vs. SVZ-targeted reprogramming

library(tidyverse)
library(biomaRt)
library(fgsea)
library(DT)
library(org.Mm.eg.db)
library(NCmisc)
library(clusterProfiler)
library(cowplot)
library(viridis)
library(pheatmap)
theme_set(theme_cowplot())

#==========
# Compare across cell types for specific pathways in whole-body and SVZ-targeted reprogramming ----
#==========
setwd("~/Dropbox/10x_expAC/4.GSEA/")

# Load data and combine
gsea_svz <- read.table("data/leadingedge/ALL_Dox.untr_fgsea_GO_2022-05-31.txt")
gsea_whole <-read.table("~/Dropbox/10x_OSKM_2/4.GSEA/data_09242020/leadingedge/ALL_OSKM-old_2Dox0.untr_fgsea_GO_2021-01-20.txt")

# Define common cell types
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "Endothelial",
           "Microglia", "Mono_Mac")

gsea_whole$celltype <- plyr::mapvalues(gsea_whole$celltype, from="Macrophage", to="Mono_Mac")
gsea_whole$celltype <- factor(gsea_whole$celltype, levels=CELLS)
gsea_svz$celltype <- factor(gsea_svz$celltype, levels=CELLS)

gsea_whole$comparison <- "whole"
gsea_svz$comparison <- "svz"

gsea_complete <- bind_rows(gsea_whole, gsea_svz)

gsea_complete$comparison <- factor(gsea_complete$comparison, levels=c("whole","svz"), ordered=T)



#============
# Scatter plots ----

widegsea <- reshape(gsea_complete, idvar = c("pathway","celltype"), timevar = "comparison", v.names=c("NES","padj","pval","log2err","size"), direction="wide")
widegsea <- subset(widegsea, select=-c(ES,leadingEdge))

# Celltype colors
CELLS2 <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
            "Oligodendrocyte", "OPC", "Endothelial",
            "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle")

data <- widegsea %>% subset(celltype %in% CELLS2)
data$celltype <- factor(data$celltype, levels=CELLS2)

newcolors <- c("#03c03c",  "#4E79A7", "#966fd6", 
               "#ffdf00", "#B6992D", "#FFBE7D", 
               "#aec6cf", "#86BCB6", "#f78d76", 
               "#db7093", "#A52A2A", "#2F4F4F")
names(newcolors) <- levels(CELLS2)


# Plots
ggplot(data %>% filter(padj.whole<0.05|padj.svz<0.05), aes(x=NES.whole, y=NES.svz)) + 
  geom_point(size=0.8, shape=16, alpha=0.6, aes(color=padj.svz)) + 
  facet_wrap(~celltype, nrow=2) +
  scale_color_viridis(limits=c(0,0.15),oob = scales::squish) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=9), axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_rect(color = "darkgray", size=1, fill=NA)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) 

ggsave("plots_wholebody_SVZtargeted/scatter_padj-svz-color_2rows_grayborder_legend.pdf", height=3.25, width=6)




#============
# Which pathways ----

# aNSCs ----
sigansc <- data %>% filter(celltype=="aNSC_NPC") %>% filter(padj.svz<0.05|padj.whole<0.05) 
sub <- gsea_complete %>% filter(celltype=="aNSC_NPC" & pathway %in% sigansc$pathway)

# Make pathway names more legible
sub$pathway <- tolower(sub$pathway)
sub$pathway <- gsub("_", " ", sub$pathway)
sub$pathway <- gsub("^...", "", sub$pathway)

# Heatmap
test <- sub %>% select(pathway, comparison, NES) %>% spread(comparison, NES)
m <- as.matrix(test[,-1])
rownames(m) <- test[,1]

p <-pheatmap(m, scale="none", cluster_cols = F, fontsize = 7, border_color = NA, treeheight_row = 15, 
             cellwidth = 20, cellheight = 5.5, cutree_rows = 4)

save_pheatmap_pdf(p, "plots_wholebody_SVZtargeted/pheatmap_whole-svz_aNSC-NPC_padj0.05either_cut4.pdf", height=4,width=4)


# neuroblast ----
signb <- data %>% filter(celltype=="Neuroblast") %>% filter(padj.svz<0.05|padj.whole<0.05) 
sub <- gsea_complete %>% filter(celltype=="Neuroblast" & pathway %in% signb$pathway)

# Make pathway names more legible
sub$pathway <- tolower(sub$pathway)
sub$pathway <- gsub("_", " ", sub$pathway)
sub$pathway <- gsub("^...", "", sub$pathway)


# Heatmap
test <- sub %>% select(pathway, comparison, NES) %>% spread(comparison, NES)
m <- as.matrix(test[,-1])
rownames(m) <- test[,1]

p <-pheatmap(m, scale="none", cluster_cols = F, fontsize = 7, border_color = NA, treeheight_row = 15, 
             cellwidth = 20, cellheight = 5.5, cutree_rows = 4)

save_pheatmap_pdf(p, "plots_wholebody_SVZtargeted/pheatmap_whole-svz_Neuroblast_padj0.05either_cut4.pdf", height=15,width=5)



#============
# Make 4 clusters of pathways for each cell type, filter to significant in at least one comparison ----

for(CELL in CELLS){

    print(CELL)
  
    # data ----
    sig <- data %>% filter(celltype==CELL) %>% filter(padj.svz<0.05|padj.whole<0.05) 
    sub <- gsea_complete %>% filter(celltype==CELL & pathway %in% sig$pathway)
    
    ## Make pathway names more legible
    sub$pathway <- tolower(sub$pathway)
    sub$pathway <- gsub("_", " ", sub$pathway)
    sub$pathway <- gsub("^...", "", sub$pathway)
    
    # pretty heatmap ----
    test <- sub %>% select(pathway, comparison, NES) %>% spread(comparison, NES)
    m <- as.matrix(test[,-1])
    rownames(m) <- test[,1]
    m <- na.omit(m)
    
    set.seed(123) #clustering will change slightly every time with random seed
    p <-pheatmap(m, scale="none", cluster_cols = F, fontsize = 7, border_color = NA, treeheight_row = 15, 
                 cellwidth = 20, cellheight = 6, kmeans_k = 4, breaks= seq(-2,2,by=4/100), main=CELL)
    p
    
    save_pheatmap_pdf(p, paste0("plots_wholebody_SVZtargeted/pheatmap_4clusters_whole-svz_",CELL,"_padj0.05either.pdf"), height=1,width=3)

}





#==========
# UpSet plots ----

library(UpSetR)

# Pull all pathways for each cell type: padj<0.05 for at least one comparison, NES>0 in both comparisons
mylist <- c()
mylist2 <- c()
mylist.next <- c()
mylist2.next <- c()
overlap1 <- c()
overlap2 <- c()
final <- c()

for (CELL in CELLS) {
 
  # significantly up in whole-body
  mylist[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% filter(padj<0.05) %>% 
    filter(comparison=="whole" & NES>0) %>% dplyr::pull(pathway)
  # up (no pval filter) in svz-targeted
  mylist.next[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% 
    filter(comparison=="svz" & NES>0) %>% dplyr::pull(pathway)
  # overlap
  overlap1[[CELL]] <- Reduce(intersect, list(mylist[[CELL]], mylist.next[[CELL]]))
  
  # significantly up in svz-targeted
  mylist2[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% filter(padj<0.05) %>% 
    filter(comparison=="svz" & NES>0) %>% dplyr::pull(pathway)
  # up (no pval filter) in whole-body
  mylist2.next[[CELL]] <-  gsea_complete %>% filter(celltype==CELL) %>% 
    filter(comparison=="whole" & NES>0) %>% dplyr::pull(pathway)
  # overlap
  overlap2[[CELL]] <- Reduce(intersect, list(mylist2[[CELL]], mylist2.next[[CELL]]))
  
  # combine overlapped lists
  final[[CELL]] <- c(overlap1[[CELL]], overlap2[[CELL]])
  
  
}


pdf("plots_wholebody_SVZtargeted/upsetplot_celltypes_pathways-NESupboth-padj0.05either_new.pdf", width = 4, height = 3)
upset(fromList(final), order.by = "freq", sets=CELLS, sets.bar.color = "darkred", mb.ratio = c(0.6,0.4),nintersects=50, text.scale = 1.2)
dev.off()


Reduce(intersect, list(final$Neuroblast,final$Microglia)) # immune/inflam

Reduce(intersect, list(final$Neuroblast,final$Microglia, final$Endothelial)) # immune




# Pull all pathways for each cell type: padj<0.05 for at least one comparison, NES<0 in both comparisons
mylist <- c()
mylist2 <- c()
mylist.next <- c()
mylist2.next <- c()
overlap1 <- c()
overlap2 <- c()
final <- c()

for (CELL in CELLS) {
  
  # significantly up in whole-body
  mylist[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% filter(padj<0.05) %>% 
    filter(comparison=="whole" & NES<0) %>% dplyr::pull(pathway)
  # up (no pval filter) in svz-targeted
  mylist.next[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% 
    filter(comparison=="svz" & NES<0) %>% dplyr::pull(pathway)
  # overlap
  overlap1[[CELL]] <- Reduce(intersect, list(mylist[[CELL]], mylist.next[[CELL]]))
  
  
  # significantly up in svz-targeted
  mylist2[[CELL]] <- gsea_complete %>% filter(celltype==CELL) %>% filter(padj<0.05) %>% 
    filter(comparison=="svz" & NES<0) %>% dplyr::pull(pathway)
  # up (no pval filter) in whole-body
  mylist2.next[[CELL]] <-  gsea_complete %>% filter(celltype==CELL) %>% 
    filter(comparison=="whole" & NES<0) %>% dplyr::pull(pathway)
  # overlap
  overlap2[[CELL]] <- Reduce(intersect, list(mylist2[[CELL]], mylist2.next[[CELL]]))
  
  # combine overlapped lists
  final[[CELL]] <- c(overlap1[[CELL]], overlap2[[CELL]])
  
  
}


pdf("plots_wholebody_SVZtargeted/upsetplot_celltypes_pathways-NESdownboth-padj0.05either_new.pdf", width = 4, height = 3)
upset(fromList(final), order.by = "freq", sets=CELLS, sets.bar.color = "darkblue", mb.ratio = c(0.6,0.4),nintersects=50, text.scale = 1.2)
dev.off()


Reduce(intersect, list(final$Neuroblast,final$Microglia, final$Oligodendrocyte)) 
#"GO_MRNA_METABOLIC_PROCESS" "GO_RNA_SPLICING"    

Reduce(intersect, list(final$Neuroblast,final$Microglia)) 
# [1] "GO_MRNA_METABOLIC_PROCESS"                    "GO_TRANSCRIPTION_REGULATOR_ACTIVITY"         
# [3] "GO_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY" "GO_RNA_SPLICING"  





#================
# Functions ----

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#================


