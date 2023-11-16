# Calculate and plot cell type proportions 

library(tidyverse)
library(Seurat)
library(ggthemes)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())

setwd("~/Dropbox/10x_OSKM_2/5.Proportions")

#========================
# Load data from cohort 1
#========================

# Read data
obj <- readRDS("../../OSKM_10x/3.Seurat/data_07012020/svz_celltypes_2020-06-29.rds")
meta <- obj[[]]

# Forget obj to free memory
rm(obj)

# Convert to tidy format
meta <- meta %>% rownames_to_column("Cell") %>% tbl_df()
meta$Age_Treatment <- plyr::mapvalues(meta$Treatment, from=c("untr","2Dox0","2Dox5"), to=c("old_untr","old_2Dox0","old_2Dox5"))
order <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")
meta$Age_Treatment <- factor(meta$Age_Treatment, levels = order, ordered = T)

# Rename
meta2 <- meta
meta2$Celltype <- meta2$Celltype.LowRes
meta2$LMO_maxID <- meta2$LMO2_maxID

#========================
# Load data from cohort 2
#========================

setwd("~/Dropbox/10x_OSKM_2/5.Proportions")

# Read data
obj <- readRDS("../3.Seurat/data_09032020/svz_subset_OSKM_2020-09-24.rds")
meta <- obj[[]]

# Forget obj to free memory
rm(obj)

# Convert to tidy format
meta <- meta %>% rownames_to_column("Cell") %>% tbl_df()
order <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")
meta$Age_Treatment <- factor(meta$Age_Treatment, levels = order, ordered = T)

#========================
# Calculate cell type proportions as a fraction of total cells per mouse
#========================

# Resolve columns and combine both experiments
colnames(meta)
colnames(meta2)

meta.1 <- meta %>% select(Cell,orig.ident,LMO_maxID,Celltype,Treatment,Age_Treatment)
meta.1$Experiment <- c("2")
meta.2 <- meta2 %>% select(Cell,orig.ident,LMO_maxID,Celltype,Treatment,Age_Treatment)
meta.2$Experiment <- c("1")

combo <- rbind(meta.1,meta.2)


# Fix factors
order <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")
combo$Age_Treatment <- factor(combo$Age_Treatment, levels = order, ordered = T)

CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Oligodendrocyte", "OPC", "Endothelial", "Microglia", "Macrophage", "Mural", "Ependymal")
combo$Celltype <- factor(combo$Celltype, levels = CELLTYPES, ordered = T)

# Calculate frequencies
Celltype.Freq <- combo %>% group_by(Age_Treatment, LMO_maxID, Experiment, Celltype) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(percent = paste0(round(100 * n/sum(n), 2), "%"))

# Calculate mean of old control mice for each cell type 
means <- Celltype.Freq %>% filter(Age_Treatment=="old_untr") %>% group_by(Age_Treatment,Celltype) %>% 
  summarise (mean_olduntr = mean(freq))
Celltype.Freq <- dplyr::inner_join(Celltype.Freq, means[,c("Celltype","mean_olduntr")], by="Celltype")

# Calculate log2 fold change over mean of old control mice
Celltype.Freq$log2FC <- log(Celltype.Freq$freq/Celltype.Freq$mean_olduntr,2)

# Subset to mice included in this study. Exclude unknown cell types and ependymal cells (see methods)
Celltype.Freq <- Celltype.Freq %>% filter(Celltype!="NA"&Celltype!="Ependymal") %>% filter(Age_Treatment!="old_2Dox5")

# Save
write.csv(Celltype.Freq, "proportions_combined-exps.csv")
Celltype.Freq <- read.csv("proportions_combined-exps.csv")

# Fix factors
order <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")
Celltype.Freq$Age_Treatment <- factor(Celltype.Freq$Age_Treatment, levels = order, ordered = T)

CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Oligodendrocyte", "OPC", "Endothelial", "Microglia", "Macrophage", "Mural", "Ependymal")
Celltype.Freq$Celltype <- factor(Celltype.Freq$Celltype, levels = CELLTYPES, ordered = T)

Celltype.Freq$Experiment <- factor(Celltype.Freq$Experiment)


#========================
# Plots
#========================

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal")
new_colors <- c("#03c03c", "#4E79A7", "#966fd6",
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#db7093", "#2F4F4F")
names(new_colors) <- CELLS


ggplot(Celltype.Freq, aes(x = Celltype, y = freq)) +
  geom_bar(aes(group = Age_Treatment, fill = Celltype), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype, y = freq, shape = Age_Treatment, color = Celltype), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = new_colors, guide='none') +
  scale_color_manual(values = new_colors, guide='none') +
  scale_shape(guide='none') +
  ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) 
ggsave("plots/combinedexps.manhattan.proportions.OSKM-no2Dox5.pdf", width = 6, height = 4)

ggplot(Celltype.Freq, aes(x = Celltype, y = log2FC)) +
  geom_bar(aes(group = Age_Treatment, fill = Celltype), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype, y = log2FC, shape = Age_Treatment, color = Celltype), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = new_colors, guide='none') +
  scale_color_manual(values = new_colors, guide='none') +
  scale_shape(guide='none') +
  ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) 
ggsave("plots/combinedexps.manhattan.proportions.log2FC.OSKM-no2Dox5.pdf", width = 6, height = 4)


# Excluding small clusters
testing <- Celltype.Freq %>% group_by(Celltype,Age_Treatment) %>% dplyr::summarize(total=sum(n))
small <- testing %>% filter(total<=30) 

ggplot(Celltype.Freq %>% filter(!Celltype %in% small$Celltype), aes(x = Celltype, y = log2FC)) +
  geom_bar(aes(group = Age_Treatment, fill = Celltype), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype, y = log2FC, shape = Age_Treatment, color = Celltype), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = new_colors, guide='none') +
  scale_color_manual(values = new_colors, guide='none') +
  scale_shape(guide='none') +
  #ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) 
ggsave("plots/combinedexps.manhattan.proportions.OSKM-no2Dox5-nosmallclusters.pdf", width = 6, height = 3.8)

#========================
# Calculate cell type proportions as a fraction of the NSC lineage per mouse
#========================

combolin <- combo %>% filter(Celltype=="Astrocyte_qNSC"|Celltype=="aNSC_NPC"|Celltype=="Neuroblast")

# Calculate frequencies
Lin.Freq.combo <- combolin %>% group_by(Age_Treatment, LMO_maxID, Experiment, Celltype) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(percent = paste0(round(100 * n/sum(n), 2), "%"))

# Save
write.csv(Lin.Freq.combo, "lineageproportions_combined-exps.csv")

doxagecolors <- c("turquoise2","firebrick", "mediumorchid")
names(doxagecolors) <- c("young_untr","old_untr","old_2Dox0")

#========================
# Boxplots

ggplot(Lin.Freq.combo %>% filter(Celltype!="Astrocyte_qNSC"), aes(x=Age_Treatment, y=freq, fill=Age_Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(shape=16) +
  facet_wrap(vars(Celltype), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = doxagecolors[1:3]) +
  theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(legend.title = element_text(size=9), 
        legend.text = element_text(size=9)) +
  ylab("Proportion of NSC lineage") +
  theme(axis.title.y = element_text(size = 9, vjust = -0.5), 
        axis.text.y = element_text(size = 8, margin = margin(0,1,0,1))) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list( c("old_2Dox0", "old_untr"), c("old_untr","young_untr")),
                     size = 2.75,tip.length = 0, step.increase=0.08) +
  scale_y_continuous(expand=expansion(mult= c(0.05, 0.08))) 
ggsave("plots/forpaper_combinedexps.lineageproportions.act-nb.y-o-2dox0.stats.pdf", width = 4.2, height = 2.4)



#-----------
# R version 4.0.2
# [1] cowplot_1.1.1      ggpubr_0.4.0       ggthemes_4.2.4     SeuratObject_4.0.2
# [5] Seurat_4.0.4       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.7       
# [9] purrr_0.3.4        readr_2.0.1        tidyr_1.1.3        tibble_3.1.4      
# [13] ggplot2_3.3.5      tidyverse_1.3.1  