# DESeq analysis of bulk RNA-seq time course, primary NSCs +/- differentiation and +/- reprogramming
# Using mapped counts (chimeric settings) as input


## **NB: Well H8 was mislabeled in sequencing run files. Should be 966 Diff 6h Dox (not untr).  
## Name of this sample is corrected in code below.

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

# Local directory. All subsequent paths are relative.
setwd("~/Dropbox/bulkRNAseq_expAL/0.Data")


# SET COLORS -----------------------------

#show_col(viridis_pal()(10))
cols <- c("#440154FF","#1F9E89FF")
names(cols) <- c("untr","Dox")

doxcolors <- c("turquoise","blue","tomato", "#612a95")
names(doxcolors) <- c("young-untr","young-Dox","old-untr","old-Dox")

timecols <- c("#000004FF", "#420A68FF", "#932667FF", "#DD513AFF", "#FCA50AFF", "#FDE725FF")
names(timecols) <- c("A06h", "A18h", "A48h", "D06h", "D18h", "D48h")

shapes <- c(1,2,5,6,15,16,17,18,19,20,7,8,9,10,3,4)


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


# =======================
# Make count matrix ----
# =======================

# Join all count files into one table ----
countFiles <- list.files("counts_chim")

counts <- data.frame(character(), numeric())
names(counts) <- c("Geneid","Length")

for (f in countFiles){
  sampleName <- strsplit(f, "_Aligned")[[1]][1]
  sampleData <- read.table(paste("counts_chim/",f,sep=""), sep="\t", dec=".", header=TRUE, stringsAsFactors=FALSE)
  names(sampleData)[7] <- sampleName
  
  counts <- full_join(counts, sampleData[,c(1,6,7)], by=c("Geneid","Length"))
}

# Correct H8
names(counts)[names(counts) == 'H8_966_D6h_untr_chim'] <- 'H8_966_D6h_Dox_chim'

str(counts)
dim(counts) #31057    67


# Calculate TPM (transcripts per kilobase million) ----
TPM <- counts
for (i in 3:ncol(counts)){
  TPM[,i] <- TPM[,i]/TPM$Length
  TPM[,i] <- TPM[,i]/sum(TPM[,i])*1000000
}
TPM$Length <- NULL
dim(TPM) #31057    66


# =======================
# Save and load ----
# =======================
#write.table(counts, "countTable.csv", sep=",", dec=".", row.names=FALSE)
#write.table(TPM, "TPMTable.csv", sep=",", dec=".", row.names=FALSE)

counts <- read.csv("countTable.csv")
TPM <- read.csv("TPMTable.csv")

# Map gene names to ensembl ids ----

counts$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                   keys=counts$Geneid, keytype = "ENSEMBL", multiVals = 'first')

TPM$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                        keys=TPM$Geneid, keytype = "ENSEMBL", multiVals = 'first')

#If no gene symbol found, put ensembl id back.
counts$SYMBOL[is.na(counts$SYMBOL)] <- counts$Geneid[is.na(counts$SYMBOL)] 
TPM$SYMBOL[is.na(TPM$SYMBOL)] <- TPM$Geneid[is.na(TPM$SYMBOL)] 

# Add gene symbols for custom transgenes
counts[31054:31057,]$SYMBOL <- c("rtTA","P2A","T2A","E2A")
TPM[31054:31057,]$SYMBOL <- c("rtTA","P2A","T2A","E2A")
tail(TPM)

# Save----
#write.table(counts, "countTable_symbols.csv", sep=",", dec=".", row.names=FALSE)
#write.table(TPM, "TPMTable_symbols.csv", sep=",", dec=".", row.names=FALSE)


# =======================
# QC: check TPM distribution across all samples----
# =======================

# Melt for plotting----
TPM_m <- melt(TPM, id.vars=c("Geneid"), variable.name="Sample", value.name="Count") #2049762       4

TPM$baseMean <- rowMeans(TPM[,2:66], na.rm=TRUE)

ggplot(TPM, aes(x=baseMean)) +
  geom_histogram(bins=100) +
  xlim(-1,10) 

TPM_m <- left_join(TPM_m, TPM[,c("Geneid", "baseMean")], by="Geneid")

# Add metadata ----
TPM_m <- separate(TPM_m, Sample, into=c("Well","MouseID","Timepoint","Treatment",NA), sep="_", remove=F)

TPM_m$Age <- plyr::mapvalues(TPM_m$MouseID, 
                             from=c("940","946","966","967","1892","1893","1894","1914","1926"),
                             to=c("old","old","old","old","young","young","young","young","young"))
TPM_m$Age <- factor(TPM_m$Age, levels=c("young","old"), ordered=T)
TPM_m$Treatment <- factor(TPM_m$Treatment, levels=c("untr","Dox"), ordered=T)
TPM_m$Timepoint <- plyr::mapvalues(TPM_m$Timepoint, 
                                   from=c("A6h","D6h"),
                                   to=c("A06h","D06h"))
TPM_m$Celltype <- factor(substr(TPM_m$Timepoint, 1, 1), levels=c("A","D"), ordered=T)
TPM_m$Hours <- factor(substr(TPM_m$Timepoint, 2, 4), levels=c("06h","18h","48h"), ordered=T)
TPM_m$Condition <- factor(paste(TPM_m$Timepoint, TPM_m$Age, TPM_m$Treatment, sep="-"))
TPM_m$Age_Treatment <- factor(paste(TPM_m$Age, TPM_m$Treatment, sep="-"), 
                              levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)
TPM_m$Timepoint_Treatment <- paste(TPM_m$Timepoint, TPM_m$Treatment, sep="-")

# Plot TPM across samples ----
ggplot(TPM_m, aes(x=Sample, y=log2(Count))) +
  geom_boxplot(outlier.size=.5) +
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="TPM")
ggsave("../1.DESeq2/plots_chim/boxplot_logTPM_sample.pdf", height = 3.5, width = 6)

ggplot(TPM_m, aes(x=Sample, y=log2(Count), col=Age_Treatment)) +
  geom_boxplot(outlier.size=.5) +
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="TPM") +
  scale_color_manual(values=doxcolors, guide='none')
ggsave("../1.DESeq2/plots_chim/boxplot_logTPM_sample_conditioncolors.pdf", height = 3.5, width = 6)



# =======================
# =======================
# DESeq2 analysis ----
# =======================
# =======================

counts <- read.csv("countTable_symbols.csv")
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


# DESeq with simple model: variable is condition (all metadata combined) ----

dds2 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata2,
                              design = ~ Condition)
dds2$Condition <- relevel(dds2$Condition, ref="A18h-old-untr")

dds2 <- DESeq(dds2)


# =======================
# Pairwise comparisons ----
# =======================

# A18h-old-Dox vs. A18h-old-untr ----
res2 <- results(dds2, contrast=c("Condition","A18h-old-Dox","A18h-old-untr"))
res2Ordered <- res2[order(res2$pvalue),]
summary(res2)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 30, 0.11%
# LFC < 0 (down)     : 39, 0.15%
# outliers [1]       : 0, 0%
# low counts [2]     : 8770, 33%
# (mean count < 5)

# store log2FC and p values in table ----
results <- as.data.frame(res2Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                   keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_A18h-old-Dox_A18h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# A06h-old-Dox vs. A18h-old-untr ----
res3 <- results(dds2, contrast=c("Condition","A06h-old-Dox","A18h-old-untr"))
res3Ordered <- res3[order(res3$pvalue),]
summary(res3)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5992, 22%
# LFC < 0 (down)     : 5820, 22%
# outliers [1]       : 0, 0%
# low counts [2]     : 5675, 21%
# (mean count < 1)

# store log2FC and p values in table ----
results <- as.data.frame(res3Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_A06h-old-Dox_A18h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# A18h-old-untr vs. A18h-young-untr ----
res <- results(dds2, contrast=c("Condition","A18h-old-untr","A18h-young-untr"))
resOrdered <- res[order(res$pvalue),]
summary(res)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1755, 6.5%
# LFC < 0 (down)     : 1203, 4.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 11349, 42%
# (mean count < 31)

# store log2FC and p values in table ----
results <- as.data.frame(resOrdered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_A18h-old-untr_A18h-young-untr.csv", sep=",", dec=".", row.names=FALSE)


# D06h-old-untr vs. D06h-young-untr ----
res4 <- results(dds2, contrast=c("Condition","D06h-old-untr","D06h-young-untr"))
res4Ordered <- res4[order(res4$pvalue),]
summary(res4)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 540, 2%
# LFC < 0 (down)     : 698, 2.6%
# outliers [1]       : 0, 0%
# low counts [2]     : 9801, 37%
# (mean count < 11)

# store log2FC and p values in table
results <- as.data.frame(res4Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D06h-old-untr_D06h-young-untr.csv", sep=",", dec=".", row.names=FALSE)


# D06h-old-Dox vs. D06h-old-untr ----
res5 <- results(dds2, contrast=c("Condition","D06h-old-Dox","D06h-old-untr"))
res5Ordered <- res5[order(res5$pvalue),]
summary(res5)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 20, 0.075%
# LFC < 0 (down)     : 11, 0.041%
# outliers [1]       : 0, 0%
# low counts [2]     : 5675, 21%
# (mean count < 1)

# store log2FC and p values in table
results <- as.data.frame(res5Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D06h-old-Dox_D06h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# D18h-old-untr vs. D18h-young-untr ----
res6 <- results(dds2, contrast=c("Condition","D18h-old-untr","D18h-young-untr"))
res6Ordered <- res6[order(res6$pvalue),]
summary(res6)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 860, 3.2%
# LFC < 0 (down)     : 1210, 4.5%
# outliers [1]       : 0, 0%
# low counts [2]     : 10317, 38%
# (mean count < 15)

# store log2FC and p values in table ----
results <- as.data.frame(res6Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D18h-old-untr_D18h-young-untr.csv", sep=",", dec=".", row.names=FALSE)


# D18h-old-Dox vs. D18h-old-untr ----
res7 <- results(dds2, contrast=c("Condition","D18h-old-Dox","D18h-old-untr"))
res7Ordered <- res7[order(res7$pvalue),]
summary(res7)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 22, 0.082%
# LFC < 0 (down)     : 8, 0.03%
# outliers [1]       : 0, 0%
# low counts [2]     : 9286, 35%
# (mean count < 7)


# store log2FC and p values in table ----
results <- as.data.frame(res7Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D18h-old-Dox_D18h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# D48h-old-untr vs. D48h-young-untr ----
res8 <- results(dds2, contrast=c("Condition","D48h-old-untr","D48h-young-untr"))
res8Ordered <- res8[order(res8$pvalue),]
summary(res8)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1552, 5.8%
# LFC < 0 (down)     : 1922, 7.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 9286, 35%
# (mean count < 7)

# store log2FC and p values in table ----
results <- as.data.frame(res8Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D48h-old-untr_D48h-young-untr.csv", sep=",", dec=".", row.names=FALSE)


# D48h-old-Dox vs. D48h-old-untr ----
res9 <- results(dds2, contrast=c("Condition","D48h-old-Dox","D48h-old-untr"))
res9Ordered <- res9[order(res9$pvalue),]
summary(res9)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 28, 0.1%
# LFC < 0 (down)     : 14, 0.052%
# outliers [1]       : 0, 0%
# low counts [2]     : 12896, 48%
# (mean count < 101)

# store log2FC and p values in table ----
results <- as.data.frame(res9Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_D48h-old-Dox_D48h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# A48h-old-Dox vs. A18h-old-untr ----
res10 <- results(dds2, contrast=c("Condition","A48h-old-Dox","A18h-old-untr"))
res10Ordered <- res10[order(res10$pvalue),]
summary(res10)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 5233, 20%
# LFC < 0 (down)     : 4838, 18%
# outliers [1]       : 0, 0%
# low counts [2]     : 6706, 25%
# (mean count < 1)

# store log2FC and p values in table ----
results <- as.data.frame(res10Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_A48h-old-Dox_A18h-old-untr.csv", sep=",", dec=".", row.names=FALSE)


# A18h-young-Dox vs. A18h-young-untr ----
res11 <- results(dds2, contrast=c("Condition","A18h-young-Dox","A18h-young-untr"))
res11Ordered <- res11[order(res11$pvalue),]
summary(res11)
# out of 26830 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 604, 2.3%
# LFC < 0 (down)     : 619, 2.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 13412, 50%
# (mean count < 149)

# store log2FC and p values in table ----
results <- as.data.frame(res11Ordered)
results$Gene <- row.names(results)
results$SYMBOL <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                         keys=results$Gene, keytype = "ENSEMBL", multiVals = 'first')
write.table(results, "../1.DESeq2/data/DESeqRes_designCondition_A18h-young-Dox_A18h-young-untr.csv", sep=",", dec=".", row.names=FALSE)



# ===================
# PCA ----
# ===================

# variance stabilizing transformation ----
vsd2 <- vst(dds2, blind=FALSE)
#write.table(assay(vsd2), "../1.DESeq2/data/VSTtable.csv",row.names = T, sep=',')

# Principal component analysis ----
# Extract for plotting
rv2 <- rowVars(assay(vsd2))
select2 <- order(rv2, decreasing = TRUE)[seq_len(min(500, length(rv2)))]
pca2 <- prcomp(t(assay(vsd2)[select2, ]))
percentVar2 <- pca2$sdev^2/sum(pca2$sdev^2)

meta2 <- as.data.frame(colData(vsd2))
d2 <- as.data.frame(pca2$x)
d2$Sample <- rownames(d2)
d2 <- inner_join(d2,meta2)

loadings <- as.data.frame(pca2$rotation) 


# Plot PCA ----

ggplot(d2, aes(x=PC1, y=PC2, color=Condition)) + 
  geom_point(size=1.5, alpha=0.7, shape=16) + 
  xlab(paste0("PC1: ",round(percentVar2[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar2[2] * 100),"% variance")) +
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8), plot.title=element_text(size=12)) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1, "lines"))

ggplot(d2, aes(x=PC1, y=PC2, color=Timepoint)) + 
  geom_point(size=1.5, alpha=0.8, shape=16) + 
  xlab(paste0("PC1 (",round(percentVar2[1] * 100),"% variance)")) +
  ylab(paste0("PC2 (",round(percentVar2[2] * 100),"% variance)")) +
  theme(axis.title=element_text(size=10),axis.text=element_text(size=10), plot.title=element_text(size=12)) +
  theme(legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.1, "lines")) +
  scale_color_manual(values=timecols)
ggsave(paste0("../1.DESeq2/plots_chim_newPCA/PCA_PC1-PC2_timecolors_new.pdf"), width = 3.5, height = 2.5)




# ===================
# Heatmaps ----
# ===================

vsd2 <- vst(dds2, blind=FALSE)

# Find 50 top variable genes ----
topVarGenes <- head(order(rowVars(assay(vsd2)),decreasing=T),50)
mat <- assay(vsd2)[ topVarGenes, ] 
mat <- mat - rowMeans(mat) # centers each gene's values across samples

# Map gene names
gene_map <- mapIds(org.Mm.eg.db, column = c('SYMBOL'), 
                    keys=rownames(mat), keytype = "ENSEMBL", multiVals = 'first')
sum(is.na(gene_map)==TRUE)

mat <- mat %>% set_rownames(gene_map)

# Heatmaps of top variable genes----
df <- as.data.frame(colData(vsd2)[,c("Timepoint","Age_Treatment")]) # Labels for heatmap
anncolors <- list(Timepoint=timecols, Age_Treatment=doxcolors)

# Reorder based on Timepoint 
df$Timepoint <- factor(df$Timepoint)
df$Age_Treatment <- factor(df$Age_Treatment, levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)
df2 <- df %>% arrange(Timepoint, Age_Treatment)

indexes <- match(rownames(df2), colnames(mat))
submat <- mat[,indexes]

# Exclude sex-linked genes
submat <- submat[!rownames(submat) %in% c("Xist","Tsix","Ddx3y","Eif2s3y","Kdm5d"),]

p <- pheatmap(submat, annotation_col=df, fontsize = 8, border_color=NA,
              show_colnames = F,
              annotation_colors = anncolors,
              cluster_cols = F, cluster_rows = T, treeheight_row = 25,
              scale = 'row')

save_pheatmap_pdf(p, "../1.DESeq2/plots_chim_newPCA/heatmap_top50-VST_all_exclude-sex_scaled_ordered.pdf", height=5.5,width=5)



#===========
# Heatmaps- top differential genes between A18h.old.Dox vs A18h.old.untr ----

#vsd2 <- vst(dds2, blind=FALSE)

# Find genes padj<0.1 in A18h.old.Dox vs A18h.old.untr----
sum(res2Ordered$padj < 0.1, na.rm=TRUE)  #69
test <- head(rownames(res2Ordered),69)

# Subset to A18h
sub18h <- coldata %>% filter(Timepoint=="A18h") %>% pull(Sample)

matsub <- assay(vsd2)[test, sub18h] 
matsub <- matsub - rowMeans(matsub) # centers each gene's values across samples

# Map gene names
gene_map <- AnnotationDbi::select(org.Mm.eg.db, columns = c('SYMBOL', 'ENSEMBL'), 
                   keys=rownames(matsub), keytype = "ENSEMBL", multiVals = 'first')

gene_map$SYMBOL[is.na(gene_map$SYMBOL)] <- gene_map$ENSEMBL[is.na(gene_map$SYMBOL)] #If no gene symbol found, put ensembl id back.

matsub <- matsub %>% set_rownames(gene_map$SYMBOL)

# New heatmap A18h----
df <- as.data.frame(colData(vsd2)[,c("Timepoint","Age_Treatment")]) # Labels for heatmap
dfsub <- df %>% subset(Timepoint=="A18h")
anncolors <- list(Timepoint=timecols, Age_Treatment=doxcolors)

# Reorder based on Age_Treatment 
dfsub$Age_Treatment <- factor(dfsub$Age_Treatment, levels=c("young-untr","young-Dox","old-untr","old-Dox"), ordered=T)
dfsub <- dfsub %>% arrange(Timepoint, Age_Treatment)

indexes <- match(rownames(dfsub), colnames(matsub))
matsub <- matsub[,indexes]

set.seed(123)
p <- pheatmap(matsub, annotation_col=dfsub, fontsize = 7, border_color=NA,
              show_colnames = F,
              annotation_colors = anncolors,
              cluster_cols = F, cluster_rows = T, treeheight_row = 10,
              scale = 'row', cutree_rows = 4, breaks= seq(-2.4,2.4,by=4.8/100))

save_pheatmap_pdf(p, "../1.DESeq2/plots_chim_newPCA/heatmap_subA18h_padj0.1-A18holdDox-A18holduntr_scaled_ordered_cut4.pdf", height=6,width=4)





#================

sessionInfo()

# R version 4.3.0 (2023-04-21)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.4.1
# 
# other attached packages:
# [1] magrittr_2.0.3              cowplot_1.1.1               limma_3.56.2                pheatmap_1.0.12            
# [5] viridis_0.6.3               viridisLite_0.4.2           biomaRt_2.56.1              clusterProfiler_4.8.1      
# [9] DOSE_3.26.1                 org.Mm.eg.db_3.17.0         AnnotationDbi_1.62.1        DESeq2_1.40.1              
# [13] SummarizedExperiment_1.30.2 Biobase_2.60.0              MatrixGenerics_1.12.2       matrixStats_1.0.0          
# [17] GenomicRanges_1.52.0        GenomeInfoDb_1.36.0         IRanges_2.34.0              S4Vectors_0.38.1           
# [21] BiocGenerics_0.46.0         reshape2_1.4.4              lubridate_1.9.2             forcats_1.0.0              
# [25] stringr_1.5.0               dplyr_1.1.2                 purrr_1.0.1                 readr_2.1.4                
# [29] tidyr_1.3.0                 tibble_3.2.1                ggplot2_3.4.2               tidyverse_2.0.0    
