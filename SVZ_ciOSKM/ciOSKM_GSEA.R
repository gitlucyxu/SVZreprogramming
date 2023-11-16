# Run basic GSEA analysis on MAST DEG data for each celltype
# SVZ-OSKM cohorts 1+2

# Credit to Stephen Turner for a very nice tutorial.
# https://stephenturner.github.io/deseq-to-fgsea/

# See end of script for package versions.
library(tidyverse)
library(biomaRt)
library(fgsea)
library(DT)
library(org.Mm.eg.db)
#library(msigdf)
library(NCmisc)
library(clusterProfiler)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())


setwd("~/Dropbox/10x_expAC/4.GSEA/")


#========================
#========================
# old vs. young ----
#========================
#========================

# Load cell type differential expression results
mast <- read.csv("../3.Seurat/data_MAST/old-v-young/de.mast_old.young_df.csv")
mast$z <- p.to.Z(mast$p_val) * sign(mast$avg_log2FC)
mast$z.adj <- p.to.Z(mast$p_val_adj) * sign(mast$avg_log2FC)

# Load pathways. Genes are as human gene symbols.
setwd("~/Dropbox/10x_OSKM_2/4.GSEA/")
pathways.go <- gmtPathways("pathways/c5.all.v7.1.symbols.gmt") #GO gene sets

# Get extra gene information, combine with DE results
gene_bitr <- bitr(mast$gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL", "ENSEMBL"), OrgDb = org.Mm.eg.db)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
      distinct() %>%
      as_tibble() %>%
      na_if("") %>% 
      na.omit()
gene_info <- inner_join(bm, gene_bitr, by = c("ensembl_gene_id" = "ENSEMBL"))
df <- inner_join(mast, gene_info, by = c("gene" = "SYMBOL"))

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")

COMPARISON <- "old.young"

gsea_all <- NULL
for (cell in CELLS) {
      
      print(cell)
      
      # Filter for celltype, reduce to just gene and DE p value statistic.
      df2 <- df %>% dplyr::filter(celltype == cell)
      print(dim(df2))
      
      df2 <- df2 %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
        na.omit() %>%
        distinct() %>%
        group_by(hsapiens_homolog_associated_gene_name) %>% 
        dplyr::summarize(stat=mean(z))
      
      # Genes ranked by Z score (not log fold change)
      ranks <- deframe(df2)
      
      # Run the Gene Set Enrichment Analysis. 
      ####### Pick the pathway set here ########
      fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
      print(head(fgseaRes))
      
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      print(head(fgseaResTidy))
      
      fgseaResTidy <- as.data.frame(fgseaResTidy)
      fgseaResTidy$`celltype` <- cell
      
      fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))
      
      # Save individual tables
      fname <- paste0("data/leadingedge/", cell, "_", COMPARISON, "_fgsea_GO_leadingedge", Sys.Date(),".txt")
      write.table(fgseaResTidy, file=fname, sep = "\t")
      gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("data/leadingedge/ALL_",COMPARISON,"_fgsea_GO_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")

# Optional: Inspect results in browser.
gsea_all %>% 
  arrange(-NES) %>% 
  DT::datatable()



#========================
#========================
# old Dox vs. old untr ----
#========================
#========================

setwd("~/Dropbox/10x_expAC/4.GSEA/")

# Load cell type differential expression results
mast <- read.csv("../3.Seurat/data_MAST/oldDox-v-olduntr/de.mast_oldDox.olduntr_df.csv")
mast$z <- p.to.Z(mast$p_val) * sign(mast$avg_log2FC)
mast$z.adj <- p.to.Z(mast$p_val_adj) * sign(mast$avg_log2FC)

# Get extra gene information, combine with DE results
gene_bitr <- bitr(mast$gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL", "ENSEMBL"), OrgDb = org.Mm.eg.db)
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()
gene_info <- inner_join(bm, gene_bitr, by = c("ensembl_gene_id" = "ENSEMBL"))
df <- inner_join(mast, gene_info, by = c("gene" = "SYMBOL"))

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")

COMPARISON <- "Dox.untr"

gsea_all <- NULL
for (cell in CELLS) {
  
  print(cell)
  
  # Filter for celltype, reduce to just gene and DE p value statistic.
  df2 <- df %>% dplyr::filter(celltype == cell)
  print(dim(df2))
  
  df2 <- df2 %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
    na.omit() %>%
    distinct() %>%
    group_by(hsapiens_homolog_associated_gene_name) %>% 
    dplyr::summarize(stat=mean(z))
  
  # Genes ranked by Z score (not log fold change)
  ranks <- deframe(df2)
  
  # Run the Gene Set Enrichment Analysis. 
  ####### Pick the pathway set here ########
  fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
  print(head(fgseaRes))
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  print(head(fgseaResTidy))
  
  fgseaResTidy <- as.data.frame(fgseaResTidy)
  fgseaResTidy$`celltype` <- cell
  
  fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))
  
  # Save individual tables
  fname <- paste0("data/leadingedge/", cell, "_", COMPARISON, "_fgsea_GO_leadingedge", Sys.Date(),".txt")
  write.table(fgseaResTidy, file=fname, sep = "\t")
  gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("data/leadingedge/ALL_",COMPARISON,"_fgsea_GO_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")

# Optional: Inspect results in browser.
gsea_all %>% 
  arrange(-NES) %>% 
  DT::datatable()


#============
# Easy plotting ----
#============

setwd("~/Dropbox/10x_expAC/4.GSEA")
gsea_all <- read.table("data/leadingedge/ALL_Dox.untr_fgsea_GO_2022-05-31.txt")
COMPARISON <- "Dox.untr"
  
CELL <- "aNSC_NPC"
data <- gsea_all %>% dplyr::filter(gsea_all$celltype == CELL) %>% dplyr::filter(padj < 0.05)

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^...", "", data$pathway)

ggplot(data, aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(CELL,"     ",COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 
ggsave(paste0("plots_forpaper/",CELL,"_",COMPARISON,"_fgsea_GO_padj0.05_dotplot_xsmall.pdf"), width = 6, height = 3.5)

# Filtered
ggplot(data %>% dplyr::filter(padj<0.03) %>% dplyr::filter(dense_rank(-abs(NES))<=20), aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(CELL,"     ",COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots_forpaper/",CELL,"_",COMPARISON,"_fgsea_GO_padj0.03_NES-top20_dotplot_small.pdf"), width = 7, height = 3)



CELL <- "Neuroblast"
data <- gsea_all %>% dplyr::filter(gsea_all$celltype == CELL) %>% dplyr::filter(padj < 0.05)

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^...", "", data$pathway)

# padj<0.01, NES top 10 and bottom 10
ggplot(data %>% dplyr::filter(padj<0.01) %>% dplyr::filter(dense_rank(NES)<=10 | dense_rank(-NES)<=10), aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(CELL,"     ",COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots_forpaper/",CELL,"_",COMPARISON,"_fgsea_GO_padj0.01_NEStop10bottom10_dotplot_small.pdf"), width = 7, height = 2)



#==========
# Compare across cell types for specific pathways in Dox-v-old and old-v-young ----
#==========
setwd("~/Dropbox/10x_expAC/4.GSEA/")

# Load data and combine
gsea_yo <-read.table("data/leadingedge/ALL_old.young_fgsea_GO_2022-05-31.txt")
gsea_dox <- read.table("data/leadingedge/ALL_Dox.untr_fgsea_GO_2022-05-31.txt")

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")

gsea_yo$celltype <- factor(gsea_yo$celltype, levels=CELLS)
gsea_dox$celltype <- factor(gsea_dox$celltype, levels=CELLS)

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_yo, gsea_dox)

gsea_complete$comparison <- factor(gsea_complete$comparison, levels=c("yo","dox"), ordered=T)


# Pick pathway and subset
PWAY <- c("GO_CELL_CELL_ADHESION", "GO_CELL_ADHESION_MOLECULE_BINDING")

data <- gsea_complete %>% dplyr::filter(pathway %in% PWAY) 

# Subset to major celltypes
CELLS2 <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle")

across <- data %>% subset(celltype %in% CELLS2)
across$celltype <- factor(across$celltype, levels=CELLS2)

across$comparison <- factor(across$comparison, levels=c("yo","dox"), ordered=T)


# Dot plots
ggplot(across , aes(x=interaction(comparison,celltype), y=pathway, size=-log(padj,10), color=NES)) + 
  geom_point(alpha = 1) + 
  scale_color_gradient2(high = "red",
                        mid = "white",
                        low = "blue") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size=10, angle=45, hjust=1), axis.title.x = element_blank(),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=8), 
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, 'cm'),
        legend.box = "horizontal")

ggsave(paste0("plots_young-old_old-Dox/celltypes_", "pathways-adhesion2","_wide", "_dotplot.pdf"), width = 8, height = 2.3)



#============
# Scatter plots ----

widegsea <- reshape(gsea_complete, idvar = c("pathway","celltype"), timevar = "comparison", v.names=c("NES","padj","pval","log2err","size"), direction="wide")
widegsea <- subset(widegsea, select=-c(ES,leadingEdge))

# Subset to major celltypes
CELLS2 <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
            "Oligodendrocyte", "OPC", "Endothelial",
            "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle")

data <- widegsea %>% subset(celltype %in% CELLS2)
data$celltype <- factor(data$celltype, levels=CELLS2)

# Celltype colors
newcolors <- c("#03c03c",  "#4E79A7", "#966fd6", 
               "#ffdf00", "#B6992D", "#FFBE7D", 
               "#aec6cf", "#86BCB6", "#f78d76", 
               "#db7093", "#A52A2A", "#2F4F4F")
names(newcolors) <- levels(CELLS2)

# New column for plotting aesthetics 
# Significant in reprogramming a = 2; otherwise a = 1
data$a <- "1"
a <- data$a
a[data$padj.dox<0.05] <- "2" 
data$a <- a

# Plot
ggplot(data %>% filter(padj.yo<0.05|padj.dox<0.05), aes(x=NES.yo, y=NES.dox, color=celltype, alpha=a)) + 
  geom_point(size=0.8, aes(alpha=a, shape=a)) + 
  facet_wrap(~celltype, nrow=3) +
  scale_color_manual(values = newcolors, guide='none') +
  scale_alpha_manual(values = c(0.3,1), guide='none') +
  scale_shape_manual(values = c(1,19), guide='none') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=9), axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_rect(color = "darkgray", size=1, fill=NA)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) 

ggsave("plots_young-old_old-Dox/scatter_celltypecolors_padj-dox-fill_nolegend_3rows_grayborder.pdf", height=6, width=7)



#==========
# UpSet plots ----

library(UpSetR)


# Pull all significant pathways for each cell type, format as a list of vectors
mylist <- c()
for (CELL in CELLS) {
  mylist[[CELL]] <- gsea_complete %>% subset(comparison=="dox") %>% subset(padj < 0.05) %>% 
    filter(celltype==CELL) %>% dplyr::pull(pathway)
}

# Limit to major cell types
SUBCELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
              "Oligodendrocyte", "Endothelial",
              "Microglia")

pdf("plots_upset/upsetplot_celltypes_pathways-dox-v-old-padj0.05_sortfreq_6celltypes.pdf", width = 6, height = 3)
upset(fromList(mylist), order.by = "freq", sets=SUBCELLS, sets.bar.color = "#612a95", mb.ratio = c(0.6,0.4)) 
dev.off()




#==============================
sessionInfo()
# 
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.6
# 
# other attached packages:
#   [1] enrichplot_1.8.1            clusterProfiler_3.16.1      org.Mm.eg.db_3.11.4        
# [4] AnnotationDbi_1.50.3        DT_0.15                     fgsea_1.14.0               
# [7] biomaRt_2.44.1              EnhancedVolcano_1.6.0       ggrepel_0.8.2              
# [10] ggpubr_0.4.0                NCmisc_1.1.6                cowplot_1.1.0              
# [13] RColorBrewer_1.1-2          viridis_0.5.1               viridisLite_0.3.0          
# [16] ggthemes_4.2.0              scales_1.1.1                sctransform_0.3            
# [19] Matrix_1.2-18               forcats_0.5.0               stringr_1.4.0              
# [22] purrr_0.3.4                 readr_1.3.1                 tidyr_1.1.2                
# [25] tibble_3.0.3                ggplot2_3.3.2               tidyverse_1.3.0            
# [28] MAST_1.14.0                 SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2
# [31] DelayedArray_0.14.1         matrixStats_0.56.0          Biobase_2.48.0             
# [34] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
# [37] S4Vectors_0.26.1            BiocGenerics_0.34.0         dplyr_1.0.2                
# [40] Seurat_3.2.1     