# Plot normalized counts for genes of interest


library(tidyverse)
library(reshape2)
library(DESeq2)
library(org.Mm.eg.db)
library(DOSE)
library(clusterProfiler)
library(biomaRt)
library(viridis)
library(pheatmap)
library(limma)
library(cowplot)
theme_set(theme_cowplot())
library(magrittr)
library(ggpubr)


setwd("~/Dropbox/bulkRNAseq_expAL/2.GSEA")

# SET COLORS -----------------------------

doxcolors <- c("turquoise","blue","tomato", "#612a95")
names(doxcolors) <- c("young-untr","young-Dox","old-untr","old-Dox")

cols <- c("gray","#1F9E89FF")
names(cols) <- c("untr","Dox")

timecols <- c("#000004FF", "#420A68FF", "#932667FF", "#DD513AFF", "#FCA50AFF", "#FDE725FF")
names(timecols) <- c("A06h", "A18h", "A48h", "D06h", "D18h", "D48h")


# Load data ----

counts <- read.csv("../0.Data/countTable_symbols.csv")
cts <- as.matrix(counts[,3:67])
row.names(cts) <- counts$Geneid

coldata <- data.frame(colnames(cts))
names(coldata) <- "Sample"

coldata <- separate(coldata, Sample, into=c("Well","MouseID","Timepoint","Treatment",NA), sep="_", remove=F)
coldata$Age <- plyr::mapvalues(coldata$MouseID, 
                               from=c("940","946","966","967","1892","1893","1894","1914","1926"),
                               to=c("old","old","old","old","young","young","young","young","young"))
coldata$Age <- factor(coldata$Age, levels=c("young","old"), ordered=T)
coldata$Treatment <- factor(coldata$Treatment, levels=c("untr","Dox"), ordered=T)
coldata$Timepoint <- plyr::mapvalues(coldata$Timepoint, 
                                     from=c("A6h","D6h"),
                                     to=c("A06h","D06h"))
coldata$Celltype <- factor(substr(coldata$Timepoint, 1, 1), levels=c("A","D"), ordered=T)
coldata$Hours <- factor(substr(coldata$Timepoint, 2, 4), levels=c("06h","18h","48h"), ordered=T)
coldata$Condition <- paste(coldata$Timepoint, coldata$Age, coldata$Treatment, sep="-")
coldata$Age_Treatment <- factor(paste(coldata$Age, coldata$Treatment, sep="-"),
                                levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)
coldata$Timepoint_Treatment <- paste(coldata$Timepoint, coldata$Treatment, sep="-")

coldata2 <- as.matrix(coldata)

dim(cts)
dim(coldata2)

dds2 <- DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata2,
                               design = ~ Condition)
dds2$Condition <- relevel(dds2$Condition, ref="A18h-old-untr")

lrtdds <- DESeq(dds2, test="LRT", reduced=~1)

gene_map2 <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=rownames(lrtdds), keytype = "ENSEMBL", multiVals = 'first')
gene_map2 <- as.data.frame(gene_map2)
gene_map2 <- rownames_to_column(gene_map2, "ENSEMBL")
colnames(gene_map2) <- c("ENSEMBL","SYMBOL")


# From top variable genes: plot counts across timepoints ----

GOI <- "Gfap"
genecounts <- NULL
genecounts <- plotCounts(dds2, gene=gene_map2 %>% filter(SYMBOL==GOI) %>% pull(ENSEMBL), intgroup=c("Condition","Age_Treatment","Timepoint","Age","Treatment","Celltype"), normalized = T, returnData = T)
# plotCounts normalizes the gene counts by default and adds a pseudocount of 0.5
genecounts$Age_Treatment <- factor(genecounts$Age_Treatment, levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)

ggplot(genecounts %>% subset(Timepoint!="A06h"&Timepoint!="A48h"&Age_Treatment!="young-Dox"), aes(x=Timepoint, y=count, color=Timepoint)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=.2,height=0), size=1.5, shape=16, alpha=0.7) +
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8),plot.title=element_text(size=10)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank()) +
  theme(strip.text = element_text(size=8)) +
  scale_color_manual(values = timecols) +
  #facet_wrap(~Age_Treatment, ncol=1) +
  labs(title=GOI, y="Normalized count") +
  scale_y_continuous(trans='log10', expand=expansion(mult = c(0.1, 0.1))) +
  theme(legend.position = 'none') +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("A18h", "D06h"), c("A18h", "D18h"), c("A18h","D48h")),
                     size = 2.5, tip.length = 0, step.increase=0.1)
ggsave(paste0("plots_genes/counts_",GOI,"_timepoint_subset_color-timepoint_log_stats.pdf"), width = 2, height = 2)

kruskal.test(genecounts %>% subset(Timepoint!="A06h"&Timepoint!="A48h"&Age_Treatment!="young-Dox") %>% pull(count), genecounts %>% subset(Timepoint!="A06h"&Timepoint!="A48h"&Age_Treatment!="young-Dox") %>% pull(Timepoint))



# OSKM: dox vs. untr----

GOI <- "Pou5f1"
genecounts <- NULL
genecounts <- plotCounts(dds2, gene=gene_map2 %>% filter(SYMBOL==GOI) %>% pull(ENSEMBL), intgroup=c("Condition","Age_Treatment","Timepoint","Age","Treatment","Celltype"), normalized = T, returnData = T)
# plotCounts normalizes the gene counts by default and adds a pseudocount of 0.5
genecounts$Age_Treatment <- factor(genecounts$Age_Treatment, levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)
genecounts$Treatment <- factor(genecounts$Treatment, levels=c("untr","Dox"), ordered=T)

ggplot(genecounts %>% subset(Condition=="A06h-old-Dox"|Condition=="A18h-old-untr"), aes(x=Treatment, y=count, color=Treatment)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=.1,jitter.height=0), size=2, shape=16, alpha=0.7) +
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8),plot.title=element_text(size=10)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank()) +
  theme(strip.text = element_text(size=8)) +
  scale_color_manual(values = cols) +
  #facet_wrap(~Celltype, nrow=1, scales='free_x') +
  labs(title=GOI, y="Normalized count") +
  scale_y_continuous(expand=expansion(mult = c(0.1, 0.1))) +
  #scale_y_continuous(trans='log10') +
  theme(legend.position = 'none') +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("untr", "Dox")),
                     size = 3, tip.length = 0, step.increase=0.08)

ggsave(paste0("../1.DESeq2/plots_differential/boxplot_counts_",GOI,"_condition_subset-A06h_treatment-colors_stats.pdf"), width = 1.25, height=2)



# Differentially expressed genes, old+OSKM vs. old, A18h ----

GOI <- "Lrfn4"
genecounts <- NULL
genecounts <- plotCounts(dds2, gene=gene_map2 %>% filter(SYMBOL==GOI) %>% pull(ENSEMBL), intgroup=c("Condition","Age_Treatment","Timepoint","Age","Treatment","Celltype"), normalized = T, returnData = T)
# plotCounts normalizes the gene counts by default and adds a pseudocount of 0.5
genecounts$Age_Treatment <- factor(genecounts$Age_Treatment, levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)

ggplot(genecounts %>% subset(Timepoint=="A18h"), aes(x=Age_Treatment, y=count, color=Age_Treatment)) +
  #geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=.1,height=0), size=2, shape=16, alpha=0.7) +
  theme(axis.title=element_text(size=10),axis.text=element_text(size=10),plot.title=element_text(size=10)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_color_manual(values = alpha(doxcolors,0.7)) +
  labs(title=GOI, y="Normalized count") +
  scale_y_continuous(trans='log10', expand=expansion(mult = c(0.1, 0.1))) +
  theme(legend.position = 'none') +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("young-untr", "young-Dox"), c("young-untr", "old-untr")),
                     size = 3.5, tip.length = 0, step.increase=0.12) +
  stat_compare_means(method = "wilcox",
                     comparisons = list(c("old-untr","old-Dox")),
                     size = 3.5, tip.length = 0, step.increase=0) +
  stat_summary(fun.data = mean_se, geom="errorbar", width = 0.2) +
  stat_summary(fun="mean", geom="crossbar", width = 0.5)

ggsave(paste0("../1.DESeq2/plots_differential/counts_",GOI,"_subsetA18h_stats_new.pdf"), width = 2, height=2)



#========
sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.4.1
# [1] ggpubr_0.6.0                magrittr_2.0.3              cowplot_1.1.1              
# [4] limma_3.56.2                pheatmap_1.0.12             viridis_0.6.3              
# [7] viridisLite_0.4.2           biomaRt_2.56.1              clusterProfiler_4.8.1      
# [10] DOSE_3.26.1                 org.Mm.eg.db_3.17.0         AnnotationDbi_1.62.1       
# [13] DESeq2_1.40.1               SummarizedExperiment_1.30.2 Biobase_2.60.0             
# [16] MatrixGenerics_1.12.2       matrixStats_1.0.0           GenomicRanges_1.52.0       
# [19] GenomeInfoDb_1.36.0         IRanges_2.34.0              S4Vectors_0.38.1           
# [22] BiocGenerics_0.46.0         reshape2_1.4.4              lubridate_1.9.2            
# [25] forcats_1.0.0               stringr_1.5.0               dplyr_1.1.2                
# [28] purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                
# [31] tibble_3.2.1                ggplot2_3.4.2               tidyverse_2.0.0   

