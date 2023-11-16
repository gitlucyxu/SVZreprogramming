# Plot and compare cell type proportions for each tissue
# Tidy-R Idiomatic Version

library(tidyverse)
library(Seurat)
library(ggthemes)
library(ggpubr)

setwd("~/Dropbox/10x_expAC/5.Proportions")

# All samples
obj <- readRDS("..//3.Seurat/data/svz_celltypes_metadata_2022-05-20.rds")
meta <- obj[[]]

# Forget obj to free memory
rm(obj)

# Convert to tidy format
meta <- meta %>% rownames_to_column("Cell") %>% tibble::as_tibble()
order <- c("young_untr","old_untr","old_Dox")
meta$age_treatment <- factor(meta$age_treatment, levels = order, ordered = T)

CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron", "Oligodendrocyte", "OPC", "Endothelial", "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
meta$Celltype.LowRes <- factor(meta$Celltype.LowRes, levels = CELLTYPES, ordered = T)

# Set colors
doxagecolors <- c("turquoise","tomato", "#612a95")
doxagecolors <- alpha(doxagecolors,0.8)
names(doxagecolors) <- unique(order)

# Modify the group_by arguments to treatment, mouse, celltype columns
Celltype.Freq <- meta %>% group_by(orig.ident, age_treatment, MULTI_classification_rescued, Celltype.LowRes) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(percent = paste0(round(100 * n/sum(n), 2), "%"))

# Add zero values for samples without any neurons
age_treatment<- c("young_untr","young_untr","young_untr","young_untr","young_untr",
                  "old_untr","old_untr","old_untr",
                  "old_Dox","old_Dox","old_Dox",
                  "old_untr","old_untr","old_untr", "young_untr")
MULTI_classification_rescued <- c("A-1108L-3.9-M-lane1-untr-TGTGATGG","B-1110L-3.9-M-lane1-untr-TCAATGGC","I-1114L-4.0-F-lane4-untr-AAGTACGC", "J-1115L-4.0-F-lane4-untr-ATTCGCAC", "K-1116L-4.0-F-lane4-untr-GAGTCGAT","C-L0111-28.1-M-lane2-untr-CTCTAGAC", "E-L0114-28.1-M-lane2-untr-AGTTGCGT", "G-L0117-28.1-F-lane2-untr-GTACCTGT","D-L0112-28.1-M-lane3-Dox-ACCAATGC", "F-L0116-28.1-F-lane3-Dox-CGAACAAG", "H-L125-27.4-M-lane3-Dox-GAAGCTTG","M-L126-27.4-M-lane5-untr-CAGTTAGG","O-L132-27.4-F-lane5-untr-AAGCAGTC","Q-L0108-28.8-M-lane5-untr-ACTCGAAG", "J-1115L-4.0-F-lane4-untr-ATTCGCAC"
)
orig.ident <- c("Lane1","Lane1","Lane4","Lane4","Lane4","Lane2","Lane2","Lane2","Lane3","Lane3","Lane3","Lane5","Lane5","Lane5","Lane4")
Celltype.LowRes <- c(rep("Neuron",14), "Ependymal")
n <- rep(0,15)
freq <- rep(0,15)
percent <- rep("0%",15)

testing<- data.frame(orig.ident,age_treatment, MULTI_classification_rescued, Celltype.LowRes,n,freq,percent)

dim(Celltype.Freq)
Celltype.Freq <- rbind(Celltype.Freq, testing)
dim(Celltype.Freq)


Celltype.Freq$Experiment <- plyr::mapvalues(Celltype.Freq$orig.ident, 
                                       from=c("Lane1","Lane2","Lane3","Lane4","Lane5","Lane6"),
                                       to=c("Exp1","Exp1","Exp1","Exp2","Exp2","Exp2"))

order <- c("young_untr","old_untr","old_Dox")
Celltype.Freq$age_treatment <- factor(Celltype.Freq$age_treatment, levels = order, ordered = T)

CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron", "Oligodendrocyte", "OPC", "Endothelial", "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
Celltype.Freq$Celltype.LowRes <- factor(Celltype.Freq$Celltype.LowRes, levels = CELLTYPES, ordered = T)


# Add log2 fold change as a column
means <- Celltype.Freq %>% filter(age_treatment=="old_untr") %>% group_by(age_treatment,Celltype.LowRes) %>%
  summarise (mean_olduntr = mean(freq))

Celltype.Freq <- dplyr::inner_join(Celltype.Freq, means[,c("Celltype.LowRes","mean_olduntr")], by="Celltype.LowRes")

Celltype.Freq$log2FC <- log(Celltype.Freq$freq/Celltype.Freq$mean_olduntr,2)


# Save and load ----
write.csv(Celltype.Freq, "proportions.csv")
setwd("~/Dropbox/10x_expAC/5.Proportions")
Celltype.Freq <- read.csv("proportions.csv")

order <- c("young_untr","old_untr","old_Dox")
Celltype.Freq$age_treatment <- factor(Celltype.Freq$age_treatment, levels = order, ordered = T)

CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron", "Oligodendrocyte", "OPC", "Endothelial", "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
Celltype.Freq$Celltype.LowRes <- factor(Celltype.Freq$Celltype.LowRes, levels = CELLTYPES, ordered = T)

# Wide format for supp tables ----
df_wide <- reshape(Celltype.Freq, idvar="MULTI_classification_rescued", timevar="Celltype.LowRes", v.names="freq", direction="wide", sep="_")
write.csv(df_wide, "wide_proportions.csv")

df_wide <- reshape(Celltype.Freq, idvar="MULTI_classification_rescued", timevar="Celltype.LowRes", v.names="log2FC", direction="wide", sep="_")
write.csv(df_wide, "wide_log2FCproportions.csv")

# Totals ----

Celltype.Freq %>% dplyr::summarize(total=sum(n))
#16291

totals <- Celltype.Freq %>% group_by(Celltype.LowRes,age_treatment) %>% dplyr::summarize(total=sum(n))

# Identify cell types with 30 or fewer total cells in a condition
small <- totals %>% filter(total<=30) 
unique(small$Celltype.LowRes)
# Neuron    OPC       Ependymal



# Plots ----
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "Mural", "Ependymal")

new_colors <- c("#03c03c", "#4E79A7", "#966fd6",
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#db7093", "#2F4F4F")
names(new_colors) <- CELLS


ggplot(Celltype.Freq %>% filter(Celltype.LowRes=="Astrocyte_qNSC"|Celltype.LowRes=="aNSC_NPC"|Celltype.LowRes=="Neuroblast"|Celltype.LowRes=="Oligodendrocyte"|Celltype.LowRes=="OPC"|Celltype.LowRes=="Endothelial"|Celltype.LowRes=="Microglia"|Celltype.LowRes=="Mono_Mac"), aes(x = Celltype.LowRes, y = log2FC)) +
  geom_bar(aes(group = age_treatment, fill = Celltype.LowRes), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype.LowRes, y = log2FC, shape = age_treatment, color = Celltype.LowRes), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = new_colors) +
  scale_color_manual(values = new_colors) +
  ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank())

ggsave("plots/manhattan.proportions.log2FC.pdf", width = 9.10, height = 4.5)


# Exclude neurons and ependymal cells -----
# DEFINE COLORS BY CELLTYPE --------------------
CELLS2 <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", 
              "Oligodendrocyte", "OPC", "Endothelial",
              "Microglia", "Mono_Mac", "T-cell", 
              "Pericyte", "VascSmoothMuscle")

cellcolors2 <- c("#03c03c",  "#4E79A7", "#966fd6", 
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#f78d76", 
                "#db7093", "#A52A2A")
# Extra colors: "#e5aa70", "#59A14F"
# neuron "#499894"
# ependymal "#2F4F4F"
cellcolors2 <- alpha(cellcolors2, 0.6)
names(cellcolors2) <- CELLS2


# Without small clusters or immune
EXCLUDE <- c("OPC","Mono_Mac","T-cell")

ggplot(Celltype.Freq %>% filter(!Celltype.LowRes %in% EXCLUDE) %>%  filter(!Celltype.LowRes %in% small$Celltype.LowRes), aes(x = Celltype.LowRes, y = log2FC)) +
  geom_bar(aes(group = age_treatment, fill = Celltype.LowRes), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype.LowRes, y = log2FC, shape = age_treatment, color = Celltype.LowRes), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = cellcolors2, guide = "none") +
  scale_color_manual(values = cellcolors2, guide = "none") +
  scale_shape(guide='none') +
  #ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank())
ggsave("plots/manhattan.proportions.log2FC.all.nosmallclusters.noimmune.pdf", width = 7, height = 4)

# OPCs and immune
ggplot(Celltype.Freq %>% filter(Celltype.LowRes %in% EXCLUDE), aes(x = Celltype.LowRes, y = log2FC)) +
  geom_bar(aes(group = age_treatment, fill = Celltype.LowRes), stat = "summary", fun = "mean",  alpha = .5, position = position_dodge2(width = 1, padding = 0.1)) +
  geom_point(aes(x = Celltype.LowRes, y = log2FC, shape = age_treatment, color = Celltype.LowRes), size = 2.5, alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = cellcolors2, guide = "none") +
  scale_color_manual(values = cellcolors2, guide = "none") +
  scale_shape(guide='none') +
  #ggtitle("Cell type proportions") +
  theme(plot.title = element_text(size=20, face = "plain")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank())

ggsave("plots/manhattan.proportions.log2FC.all.OPC-immune.pdf", width = 3.5, height = 4)



#========
# aNSC/NPC and neuroblasts, as a proportion of NSC lineage ----
#========

# Modify the group_by arguments to treatment, mouse, celltype columns
Lin.Freq1 <- meta %>% filter(Celltype.LowRes=="Astrocyte_qNSC"|Celltype.LowRes=="aNSC_NPC"|Celltype.LowRes=="Neuroblast") %>%
  group_by(age_treatment, orig.ident, MULTI_classification_rescued, Celltype.LowRes) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  mutate(percent = paste0(round(100 * n/sum(n), 2), "%")) 

Lin.Freq1$Experiment <- plyr::mapvalues(Lin.Freq1$orig.ident, 
                                       from=c("Lane1","Lane2","Lane3","Lane4","Lane5","Lane6"),
                                       to=c("Exp1","Exp1","Exp1","Exp2","Exp2","Exp2"))

write.csv(Lin.Freq1, "lineageproportions.csv")

# Wide format for supp tables ----
Lin.Freq1 <- as.data.frame(Lin.Freq1)
df_wide <-  reshape(Lin.Freq1, idvar="MULTI_classification_rescued", timevar="Celltype.LowRes", v.names="freq", direction="wide", sep="_")
write.csv(df_wide, "wide_lineageproportions.csv")


# Plots
ggplot(Lin.Freq1 %>% filter(Celltype.LowRes!="Astrocyte_qNSC"), aes(x=age_treatment, y=freq, fill=age_treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(show.legend = F, alpha=0.5, position=position_jitter(width=0.05)) +
  facet_wrap(vars(Celltype.LowRes), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = doxagecolors) +
  theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(strip.text.x = element_text(size = 11)) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=9)) +
  ylab("Proportion of NSC lineage") +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("old_Dox", "old_untr"), c("old_untr","young_untr")),
                     size = 2.8,tip.length = 0, step.increase=0.07)  + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
ggsave("plots/lineage.proportions.dots.stats.pdf", width = 6.5, height = 2.5)









