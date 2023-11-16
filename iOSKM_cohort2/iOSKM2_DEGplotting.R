# Manhattan-like plots highlighting DEGs by Cell type using MAST DE values.
# Cell type proportions by sample.
# Lucy Xu 9/9/2020

rm(list = ls())
library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc) 
# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/10x_OSKM_2/3.Seurat")


#===========================================================================================
# Reprogramming ----
#===========================================================================================
# Old+OSKM vs. old control

df_2Dox0_untr <- read_csv("data_09032020/de.mast.OSKM-old.2Dox0.untr_df.csv")

# Add z-score based on two sided null hypothesis.
df_2Dox0_untr$z <- p.to.Z(df_2Dox0_untr$p_val) * sign(df_2Dox0_untr$avg_logFC)

# Randomize rows to reduce overplotting issues.
df_2Dox0_untr <- df_2Dox0_untr[sample(nrow(df_2Dox0_untr)), ]

# Adjust factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal", "27", "28","29","30")
df_2Dox0_untr$celltype_factor <- factor(df_2Dox0_untr$celltype,  levels=CELLS, ordered=T)

# Subset to known cell types and larger clusters
df_2Dox0_untr <- df_2Dox0_untr %>% filter(!celltype %in% c("27","28","29"))

# Color based on padj<0.2
col <- newcolors[df_2Dox0_untr$celltype_factor]
col[df_2Dox0_untr$p_adj_fdr > 0.2] <- "#D3D3D3" #grey
df_2Dox0_untr$col <- as.factor(col)

# Plotting
p <- ggplot(df_2Dox0_untr, aes(x = celltype_factor, y = z, color = col)) +
  geom_jitter(shape=16, width = 0.40, alpha = .55, size = 1) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("2Dox0 vs. untr") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df_2Dox0_untr$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) 
  
png("plots.manhattan_09032020/manhattan.OSKM-old.2Dox0.untr.padjfdr.0.2.new.png", width=5.5, height=5, units = 'in', res=600)
p
dev.off()


#===========================================================================================
# Aging ----
#===========================================================================================
# Old control vs. young control

df_old_young <- read_csv("data_09032020/de.mast.OSKM-untr.old.young_df.csv")

# Add z-score based on two sided null hypothesis.
df_old_young$z <- p.to.Z(df_old_young$p_val) * sign(df_old_young$avg_logFC)
df_old_young$z.adj <- p.to.Z(df_old_young$p_val_adj) * sign(df_old_young$avg_logFC)

# Randomize rows to reduce overplotting issues.
df_old_young <- df_old_young[sample(nrow(df_old_young)), ]

# Adjust factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal", "27", "28","29","30")
df_old_young$celltype_factor <- factor(df_old_young$celltype,  levels=CELLS, ordered=T)

# Subset to known cell types and larger clusters
df_old_young <- df_old_young %>% filter(!celltype %in% c("27","28","29","30","OPC","Ependymal"))

# Color based on padj<0.2
col <- newcolors[df_old_young$celltype_factor]
col[df_old_young$p_adj_fdr > 0.2] <- "#D3D3D3" # grey
df_old_young$col <- as.factor(col)

# Plotting
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

png("plots.manhattan_09032020/manhattan.OSKM-untr.young.old.padjfdr.0.2.new.png", width=5.5, height=5, units = 'in', res=600)
p
dev.off()

#======
sessionInfo()

#R version 4.0.2 (2020-06-22)
# [1] NCmisc_1.1.6       scales_1.1.1       viridis_0.6.1      viridisLite_0.4.0 
# [5] RColorBrewer_1.1-2 cowplot_1.1.1      forcats_0.5.1      stringr_1.4.0     
# [9] dplyr_1.0.7        purrr_0.3.4        readr_2.0.1        tidyr_1.1.3       
# [13] tibble_3.1.4       ggplot2_3.3.5      tidyverse_1.3.1   

