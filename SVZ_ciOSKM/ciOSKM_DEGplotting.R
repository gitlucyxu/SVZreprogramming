# Manhattan-like plots highlighting DEGs by cell type using MAST DE values.
# SVZ-OSKM cohorts 1+2

rm(list = ls())
library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc) 
theme_set(theme_cowplot())

# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/10x_expAC/3.Seurat")


#==================================================================================================
# old untr vs. young untr ----
#==================================================================================================

df_old_young <- read_csv("data_MAST/old-v-young/de.mast_old.young_df.csv")

# Add z-score based on two sided null hypothesis.
df_old_young$z <- p.to.Z(df_old_young$p_val) * sign(df_old_young$avg_log2FC)
df_old_young$z.adj <- p.to.Z(df_old_young$p_val_adj) * sign(df_old_young$avg_log2FC)

# Randomize rows to reduce overplotting issues.
df_old_young <- df_old_young[sample(nrow(df_old_young)), ]

# Adjust factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
df_old_young$celltype_factor <- factor(df_old_young$celltype,  levels=CELLS, ordered=T)

# Adjust colors
newcolors <- c("#03c03c",  "#4E79A7", "#966fd6", 
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#f78d76", 
                "#db7093", "#A52A2A", "#2F4F4F")
names(newcolors) <- levels(CELLS)

#==================================================================================================
# P_ADJ_FDR < 0.2

df_old_young <- df_old_young %>% filter(celltype!="Ependymal")

# Make custom color column to facilitate grey coloring by threshold.
col <- newcolors[df_old_young$celltype_factor]
col[df_old_young$p_adj_fdr > 0.2] <- "#D3D3D3" # grey
df_old_young$col <- as.factor(col)

p <- ggplot(df_old_young, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(shape=16, width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("old vs. young") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df_old_young$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) 
p

png("plots_manhattan/manhattan_old.young_padjfdr.0.2.new.png", width=6, height=5, units = 'in', res=600)
p
dev.off()



#==================================================================================================
# old Dox vs. old untr ----
#==================================================================================================

df_Dox_untr <- read_csv("data_MAST/oldDox-v-olduntr/de.mast_oldDox.olduntr_df.csv")

# Add z-score based on two sided null hypothesis.
df_Dox_untr$z <- p.to.Z(df_Dox_untr$p_val) * sign(df_Dox_untr$avg_log2FC)
df_Dox_untr$z.adj <- p.to.Z(df_Dox_untr$p_val_adj) * sign(df_Dox_untr$avg_log2FC)

# Randomize rows to reduce overplotting issues.
df_Dox_untr <- df_Dox_untr[sample(nrow(df_Dox_untr)), ]

# Adjust factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
df_Dox_untr$celltype_factor <- factor(df_Dox_untr$celltype,  levels=CELLS, ordered=T)

# Adjust colors
newcolors <- c("#03c03c",  "#4E79A7", "#966fd6", 
               "#ffdf00", "#B6992D", "#FFBE7D", 
               "#aec6cf", "#86BCB6", "#f78d76", 
               "#db7093", "#A52A2A", "#2F4F4F")
names(newcolors) <- levels(CELLS)

#==================================================================================================
# P_ADJ_FDR < 0.2

df_Dox_untr <- df_Dox_untr %>% filter(celltype!="Ependymal")

# Make custom color column to facilitate grey coloring by threshold.
col <- newcolors[df_Dox_untr$celltype_factor]
col[df_Dox_untr$p_adj_fdr > 0.2] <- "#D3D3D3" # grey
df_Dox_untr$col <- as.factor(col)

p <- ggplot(df_Dox_untr, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(shape=16, width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("old Dox vs. old control") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df_Dox_untr$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) 
p

png("plots_manhattan/manhattan_oldDox.olduntr_padjfdr.0.2.new.png", width=6, height=5, units = 'in', res=600)
p
dev.off()





#-------------------------------------------------------------------------

sessionInfo()

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# [1] NCmisc_1.1.5                cowplot_1.1.1               RColorBrewer_1.1-2         
# [4] viridis_0.6.1               viridisLite_0.4.0           ggthemes_4.2.4             
# [7] scales_1.1.1                sctransform_0.3.2           Matrix_1.3-4               
# [10] forcats_0.5.1               stringr_1.4.0               purrr_0.3.4                
# [13] readr_2.0.1                 tidyr_1.1.3                 tibble_3.1.4               
# [16] ggplot2_3.3.5               tidyverse_1.3.1             MAST_1.14.0                
# [19] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2 DelayedArray_0.14.1        
# [22] matrixStats_0.60.1          Biobase_2.48.0              GenomicRanges_1.40.0       
# [25] GenomeInfoDb_1.24.2         IRanges_2.22.2              S4Vectors_0.26.1           
# [28] BiocGenerics_0.34.0         dplyr_1.0.7                 SeuratObject_4.0.2         
# [31] Seurat_4.0.4            