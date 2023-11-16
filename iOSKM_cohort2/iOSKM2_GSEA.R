# Run basic GSEA analysis on MAST DEG data for each celltype

# Credit to Stephen Turner for a very nice tutorial.
# https://stephenturner.github.io/deseq-to-fgsea/

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
library(UpSetR)
theme_set(theme_cowplot())


setwd("~/Dropbox/10x_OSKM_2/4.GSEA/") # Local

#========================
# GSEA, old+OSKM vs. old control ----
#========================

# Load cell type differential expression results
mast <- read.csv("../3.Seurat/data_09032020/de.mast.OSKM-old.2Dox0.untr_df.csv")
mast$z <- p.to.Z(mast$p_val) * sign(mast$avg_logFC)

# Load GO pathways. Genes are as human gene symbols.
pathways.go <- gmtPathways("pathways/c5.all.v7.1.symbols.gmt") #GO gene sets from MSigDB

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
              "Oligodendrocyte", "Endothelial",
              "Microglia", "Macrophage", "Mural")

COMPARISON <- "OSKM-old_2Dox0.untr"

# Run GSEA
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
        dplyr::select(-leadingEdge, -ES) %>% 
        as_tibble() %>%
        arrange(desc(NES))
      print(head(fgseaResTidy))
      
      fgseaResTidy <- as.data.frame(fgseaResTidy)
      fgseaResTidy$`celltype` <- cell
      
      # Save individual tables
      fname <- paste0("data_09242020/", cell, "_", COMPARISON, "_fgsea_GO", Sys.Date(),".txt")
      #write.table(fgseaResTidy, file=fname, sep = "\t")
      gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("data_09242020/ALL_",COMPARISON,"_fgsea_GO_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")


#========================
# Dot plots ----
#========================

setwd("~/Dropbox/10x_OSKM_2/4.GSEA")
gsea_all <- read.table("data_09242020/ALL_OSKM-old_2Dox0.untr_fgsea_GO_2021-01-20.txt")
COMPARISON <- "OSKM-old_2Dox0.untr"

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "Endothelial",
           "Microglia", "Macrophage", "Mural")

for (CELL in CELLS) {
  
  data <- gsea_all %>% dplyr::filter(gsea_all$celltype == CELL) %>% dplyr::filter(padj < 0.05)

  ggplot(data, aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
    geom_point(alpha = 1) + 
    scale_color_viridis(limits=c(0,0.055),oob = scales::squish) +
    theme_minimal() +
    labs(x="Normalized Enrichment Score", 
         subtitle=paste0(CELL,"     ",COMPARISON)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
          axis.title.y=element_blank()) +
    theme(legend.title = element_text(size=10), 
          legend.text = element_text(size=10)) 
  
  ggsave(paste0("plots_forpaper/",CELL,"_",COMPARISON,"_fgsea_GO_padj0.05_dotplot_small.pdf"), width = 10, height = 5)
}


#========================
# GSEA, old+OSKM vs. old control ----
#========================

# Load cell type differential expression results
mast <- read.csv("../3.Seurat/data_09032020/de.mast.OSKM-untr.old.young_df.csv")
mast$z <- p.to.Z(mast$p_val) * sign(mast$avg_logFC)

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
           "Oligodendrocyte", "OPC","Endothelial",
           "Microglia", "Macrophage", "Mural","Ependymal")

COMPARISON <- "OSKM.old.young"

# Run GSEA
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
    summarize(stat=mean(z))
  
  # Genes ranked by Z score (not log fold change)
  ranks <- deframe(df2)
  
  # Run the Gene Set Enrichment Analysis. 
  ## Pick the pathway set here. 
  fgseaRes <- fgsea(pathways=pathways.go, stats=ranks)
  print(head(fgseaRes))
  
  fgseaResTidy <- fgseaRes %>%
    dplyr::select(-leadingEdge, -ES) %>% 
    as_tibble() %>%
    arrange(desc(NES))
  print(head(fgseaResTidy))
  
  fgseaResTidy <- as.data.frame(fgseaResTidy)
  fgseaResTidy$`celltype` <- cell
  
  # Save individual tables
  fname <- paste0("data_09242020/", cell, "_", COMPARISON, "_fgsea_GO", Sys.Date(),".txt")
  write.table(fgseaResTidy, file=fname, sep = "\t")
  gsea_all <- rbind(gsea_all, fgseaResTidy)
}

# Save complete dataframe
fname <- paste0("data_09242020/ALL_",COMPARISON,"_fgsea_GO_", Sys.Date(),".txt")
write.table(gsea_all, file=fname, sep = "\t")



#==========
# Compare young vs. old vs. OSKM ----
#==========

setwd("~/Dropbox/10x_OSKM_2/4.GSEA")

# Load data and combine
gsea_yo <-read.table("data_09242020/ALL_OSKM.old.young_fgsea_GO_2020-09-29.txt")
gsea_dox <- read.table("data_09242020/leadingedge/ALL_OSKM-old_2Dox0.untr_fgsea_GO_2021-01-20.txt")

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "Endothelial",
           "Microglia", "Macrophage", "Mural")

gsea_yo$celltype <- factor(gsea_yo$celltype, levels=CELLS)
gsea_dox$celltype <- factor(gsea_dox$celltype, levels=CELLS)

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_yo, gsea_dox)

gsea_complete$comparison <- factor(gsea_complete$comparison, levels=c("yo","dox"), ordered=T)


#----------
# Scatter plots ----
#----------

widegsea <- reshape(gsea_complete, idvar = c("pathway","celltype"), timevar = "comparison", v.names=c("NES","padj","pval","log2err","size"), direction="wide")

# Subset to known celltypes
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "Endothelial",
           "Microglia", "Macrophage", "Mural")

data <- widegsea %>% subset(celltype %in% CELLS)
data$celltype <- factor(data$celltype, levels=CELLS)

# Make custom columns to facilitate coloring and/or other aesthetics by threshold.
newcolors <- c("#03c03c", "#4E79A7", "#966fd6",
               "#ffdf00",  "#FFBE7D", "#aec6cf", "#86BCB6", 
               "#e5aa70")
names(newcolors) <- levels(CELLS)

col <- newcolors[data$celltype]
data$col <- as.factor(col)

# Set a=1 for pathways that are significant in reprogramming, a=0.6 if not
data$a <- "1"
a <- data$a
a[data$padj.dox>0.05] <- "0.6" 
data$a <- a

# Plot, filtered to pathways that are significant in at least one comparison
ggplot(data %>% filter(padj.yo<0.05|padj.dox<0.05), aes(x=NES.yo, y=NES.dox, color=celltype, alpha=a)) + 
  geom_point(size=0.8, aes(alpha=a, shape=a)) + 
  facet_wrap(~celltype, nrow=3) +
  scale_color_manual(values = newcolors, guide='none') +
  scale_alpha_manual(values = c(0.5,1), guide='none') +
  scale_shape_manual(values = c(1,19), guide='none') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=9), axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_rect(color = "darkgray", size=1, fill=NA)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) 

ggsave("plots_young-old_old-2dox0/scatter_celltypecolors_padj-dox-fill_nolegend_3rows_grayborder.pdf", height=5, width=5)



#----------
# Compare specific pathways across cell types ----
#----------

# Pick pathways and subset
PWAY <- c("GO_RNA_PROCESSING","GO_RNA_BINDING","GO_RIBOSOME","GO_CYTOSOLIC_RIBOSOME")

data <- gsea_complete %>% dplyr::filter(pathway %in% PWAY) 

# Subset to major celltypes
across <- data %>% subset(celltype %in% CELLS)
across$celltype <- factor(across$celltype, levels=CELLS)

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
        legend.key.size = unit(0.2, 'cm'),
        legend.box = "horizontal")

ggsave(paste0("plots_young-old_old-Dox/celltypes_", "pathways-ribo","_wide", "_dotplot.pdf"), width = 8.5, height = 2.5)


#----------
# Compare pathways changing in aging vs. reprogramming within a cell type ----
#----------

# Pick cell type and subset
CELL <- "Endothelial"

data <- gsea_complete %>% dplyr::filter(celltype == CELL) 

# Subset to padj < 0.05
sub <- data %>% filter(padj < 0.05)
sig <- data %>% subset(pathway %in% sub$pathway)

# Reorder pathways for plotting based on NES of young-v-old
vals <- sig %>% subset(comparison=="yo") %>% arrange(desc(NES)) %>% 
  mutate(rank = rank(NES))
vals <- vals[,c("pathway","rank")]  

sig <- merge(sig, vals, by="pathway")  

## Make pathway names more legible
sig$pathway <- tolower(sig$pathway)
sig$pathway <- gsub("_", " ", sig$pathway)
sig$pathway <- gsub("^...", "", sig$pathway)

# Dot plots
ggplot(sig, aes(x=comparison, y=reorder(pathway, rank), size=-log(padj,10), color=NES)) + 
  geom_point(alpha = 1) + 
  scale_color_gradient2(high = "red",
                        mid = "white",
                        low = "blue") +
  theme_minimal() +
  labs(subtitle=CELL) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size=10), 
        axis.title =element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots_young-old_old-2dox0/", CELL, "_fgsea_GO_padj0.05_dotplot_new.pdf"), width = 6.5, height = 10)




#==========
# UpSet plots ----
#==========

# Pull all significant pathways for each cell type, format as a list of vectors
mylist <- c()
for (CELL in CELLS) {
  mylist[[CELL]] <- gsea_complete %>% subset(comparison=="dox") %>% subset(padj < 0.05) %>% 
    filter(celltype==CELL) %>% dplyr::pull(pathway)
}

# Upset plots

pdf("plots_upset/upsetplot_celltypes_pathways-dox-v-old-padj0.05_sortfreq.pdf", width = 10, height = 4)
upset(fromList(mylist), order.by = "freq", sets=CELLS, sets.bar.color = "#612a95", mb.ratio = c(0.6,0.4),nintersects=1000)
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
# [40] Seurat_3.2.1                UpSetR_1.3.0           