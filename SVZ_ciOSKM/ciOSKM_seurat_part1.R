# SVZ-targeted partial reprogramming (c+iOSKM))
# LMO + SVZ Seurat Processing
# 2 cohorts across 6 total 10x lanes

rm(list = ls())
library(Seurat)
library(deMULTIplex)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
library(cowplot)
theme_set(theme_cowplot())
sessionInfo()


# All subsequent paths will be relative.
setwd("~/Dropbox/10x_expAC/3.Seurat")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================
# Mapped to custom genome (mm10, ROSA26 locus: NeoR-polyAs, rtTA-IRES-EGFP-pA with updated rtTA sequence, Cre, 2A sequences)
# Cell ranger for GEX1/2/3/4/5/6
# Using LMO counts from LMO1/2/3/4/5/6
#=======================================================================
#=======================================================================
# 10x Lane 1 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data1 <- Read10X("../1.CellRanger/output1/filtered_feature_bc_matrix") #31059 x 3219, genes x cells

# Read in LMO count matrix
lmo.data1 <- Read10X("../2.LMO/output1/umi_count/", gene.column=1) #3 x 3219, samples and unmapped x cells
#lmo.data1 <- lmo.data1[1:2,] #2 x 3227, samples x cells

# Remove any cells without barcodes and barcodes without cells
colnames(svz.data1) <- sub("-1","",colnames(svz.data1))
barcode.intersect <- intersect(colnames(svz.data1), colnames(lmo.data1)) #3219
lmo.data1 <- lmo.data1[ , barcode.intersect] #3  3219
svz.data1 <- svz.data1[ , barcode.intersect] #31059  3219

# Make Seurat Object
svz1 <- CreateSeuratObject(counts = svz.data1, project = "Lane1") #31059 features across 3219 samples within 1 assay
svz1 <- NormalizeData(svz1)
svz1 <- FindVariableFeatures(svz1)
svz1 <- ScaleData(svz1)


#======================================================================================
# QC Metrics
#======================================================================================

svz1[["percent.mt"]] <- PercentageFeatureSet(svz1, pattern = "^mt-")

# cut off of less than 10% mito; more than 500 features.
svz1 <- subset(svz1, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 3012 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane1.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 18 PCs for tSNE visualization and graph-based clustering
svz1 <- RunPCA(svz1, verbose = FALSE)
ElbowPlot(svz1, ndims = 50)
ggsave("plots1/elbow.lane1.pdf", width=7.3, height = 4.3)

svz1 <- FindNeighbors(svz1, dims = 1:18)
svz1 <- FindClusters(svz1, resolution = 0.15)
svz1 <- RunUMAP(svz1, dims = 1:18)
DimPlot(svz1)
ggsave("plots1/umap.initial.fil.lane1.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

#MULTIseq Lane 1 ----

barcode.intersect <- intersect(colnames(svz1), colnames(lmo.data1)) #3012
lmo.data1 <- lmo.data1[1:2 , barcode.intersect] #3 3012
svz1[["LMO"]] <- CreateAssayObject(counts = lmo.data1)

svz1 <- NormalizeData(object = svz1, assay = "LMO", normalization.method = "CLR", margin=2)

multi <- MULTIseqDemux(svz1, assay = "LMO", autoThresh = T, maxiter = 5, 
                       qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi$MULTI_ID)
# A-1108L-3.9-M-lane1-untr-TGTGATGG B-1110L-3.9-M-lane1-untr-TCAATGGC 
# 1040                              1018 
# Negative 
# 954 

# Semi-supervised negative barcode reclassification ----
df <- t(as.matrix(multi@assays$LMO@data))
bar.table.full <- df[,1:2]
df2 <- as.data.frame(multi$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 50) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]


multi$MULTI_classification_rescued <- final.calls.rescued
table(multi$MULTI_classification_rescued)
# A-1108L-3.9-M-lane1-untr-TGTGATGG B-1110L-3.9-M-lane1-untr-TCAATGGC 
# 1422                              1428 
# Negative 
# 162 
table(multi$MULTI_classification)

Idents(multi) <- "MULTI_classification_rescued"

multi$MULTI_classification.global <- plyr::mapvalues(multi$MULTI_classification_rescued,
                                                     from = c("A-1108L-3.9-M-lane1-untr-TGTGATGG",
                                                              "B-1110L-3.9-M-lane1-untr-TCAATGGC",
                                                              "Negative","Doublet"),
                                                     to = c("Singlet", "Singlet",
                                                            "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi$MULTI_classification.global <- factor(multi$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi) <- "MULTI_classification.global"
VlnPlot(multi, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
ggsave("plots1/MULTI.CLR.2/violin.rnaCounts.globalMULTI.rescued.lane1.pdf", width = 6, height = 5)

VlnPlot(multi, features = "nFeature_RNA", pt.size = 0.1, log = F)
ggsave("plots1/MULTI.CLR.2/violin.featureCounts.globalMULTI.rescued.lane1.pdf", width = 6, height = 5)

VlnPlot(multi, features = "nFeature_RNA", pt.size = 0.1, log = T)
ggsave("plots1/MULTI.CLR.2/violin.featureCounts-log.globalMULTI.rescued.lane1.pdf", width = 6, height = 5)


#===================================
# Intermediate Saving/Loading Point ----
saveRDS(svz1, "data/svz1.05192022.rds")
saveRDS(multi, "data/svz1.MULTI.CLR.2.rescued.05192022.rds")



#=======================================================================
#=======================================================================
# 10x Lane 2 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data2 <- Read10X("../1.CellRanger/output2/filtered_feature_bc_matrix") #31059  3009, genes x cells

lmo.data2 <- Read10X("../2.LMO/output2/umi_count/", gene.column=1) #4 3009, samples and unmapped x cells

colnames(svz.data2) <- sub("-1","",colnames(svz.data2))
barcode.intersect <- intersect(colnames(svz.data2), colnames(lmo.data2)) #3009
lmo.data2 <- lmo.data2[ , barcode.intersect] 
svz.data2 <- svz.data2[ , barcode.intersect] 

svz2 <- CreateSeuratObject(counts = svz.data2, project = "Lane2") #31059 features across 3009 samples within 1 assay
svz2 <- NormalizeData(svz2)
svz2 <- FindVariableFeatures(svz2)
svz2 <- ScaleData(svz2)


#======================================================================================
# QC Metrics
#======================================================================================

svz2[["percent.mt"]] <- PercentageFeatureSet(svz2, pattern = "^mt-")

# cut off of less than 10% mito; more than 500 features.
svz2 <- subset(svz2, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 2792 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane2.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select PCs for tSNE visualization and graph-based clustering
svz2 <- RunPCA(svz2, verbose = FALSE)
ElbowPlot(svz2, ndims = 50)
ggsave("plots1/elbow.lane2.pdf", width=7.3, height = 4.3)

svz2 <- FindNeighbors(svz2, dims = 1:19)
svz2 <- FindClusters(svz2, resolution = 0.15)
svz2 <- RunUMAP(svz2, dims = 1:19)
DimPlot(svz2)
ggsave("plots1/umap.initial.fil.lane2.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz2), colnames(lmo.data2)) #2792
lmo.data2 <- lmo.data2[1:3 , barcode.intersect] 
svz2[["LMO"]] <- CreateAssayObject(counts = lmo.data2)

#MULTIseq
svz2 <- NormalizeData(object = svz2, assay = "LMO", normalization.method = "CLR", margin=2)

multi2 <- MULTIseqDemux(svz2, assay = "LMO", 
                        autoThresh = T, maxiter = 4, 
                       qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi2$MULTI_ID)
# C-L0111-28.1-M-lane2-untr-CTCTAGAC  Doublet 
# 80                                  45 
# E-L0114-28.1-M-lane2-untr-AGTTGCGT G-L0117-28.1-F-lane2-untr-GTACCTGT 
# 188                                121 
# Negative 
# 2358 

# Semi-supervised negative barcode reclassification
df <- t(as.matrix(multi2@assays$LMO@data))
bar.table.full <- df[,1:3]
df2 <- as.data.frame(multi2$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi2$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 7) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]


multi2$MULTI_classification_rescued <- final.calls.rescued
table(multi2$MULTI_classification_rescued)
# C-L0111-28.1-M-lane2-untr-CTCTAGAC  Doublet 
# 712                                 45 
# E-L0114-28.1-M-lane2-untr-AGTTGCGT G-L0117-28.1-F-lane2-untr-GTACCTGT 
# 209                                410 
# Negative 
# 1416 

Idents(multi2) <- "MULTI_classification_rescued"

multi2$MULTI_classification.global <- plyr::mapvalues(multi2$MULTI_classification_rescued,
                                                from = c("C-L0111-28.1-M-lane2-untr-CTCTAGAC",
                                                         "E-L0114-28.1-M-lane2-untr-AGTTGCGT",
                                                        "G-L0117-28.1-F-lane2-untr-GTACCTGT",
                                                              "Negative","Doublet"),
                                                     to = c("Singlet", "Singlet", "Singlet",
                                                            "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi2$MULTI_classification.global <- factor(multi2$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi2) <- "MULTI_classification.global"
VlnPlot(multi2, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

VlnPlot(multi2, features = "nFeature_RNA", pt.size = 0.1, log = F)

VlnPlot(multi2, features = "nFeature_RNA", pt.size = 0.1, log = T)

#===================================
# Intermediate Saving/Loading Point
saveRDS(svz2, "data/svz2.05192022.rds")
saveRDS(multi2, "data/svz2.MULTI.CLR.2.rescued.05192022.rds")



#=======================================================================
#=======================================================================
# 10x Lane 3 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data3 <- Read10X("../1.CellRanger/output3/filtered_feature_bc_matrix") #31059  1581, genes x cells

lmo.data3 <- Read10X("../2.LMO/output3/umi_count/", gene.column=1) # 4 1581, samples and unmapped x cells
lmo.data3 <- lmo.data3[1:3,] #3 1581, samples x cells

colnames(svz.data3) <- sub("-1","",colnames(svz.data3))
barcode.intersect <- intersect(colnames(svz.data3), colnames(lmo.data3)) #1581
lmo.data3 <- lmo.data3[ , barcode.intersect] 
svz.data3 <- svz.data3[ , barcode.intersect] 

svz3 <- CreateSeuratObject(counts = svz.data3, project = "Lane3") #31059 features across 1581 samples within 1 assay
svz3 <- NormalizeData(svz3)
svz3 <- FindVariableFeatures(svz3)
svz3 <- ScaleData(svz3)


#======================================================================================
# QC Metrics
#======================================================================================

svz3[["percent.mt"]] <- PercentageFeatureSet(svz3, pattern = "^mt-")

# cut off of less than 10% mito; more than 500 features.
svz3 <- subset(svz3, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 1462 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane3.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select PCs for tSNE visualization and graph-based clustering
svz3 <- RunPCA(svz3, verbose = FALSE)
ElbowPlot(svz3, ndims = 50)
ggsave("plots1/elbow.lane3.pdf", width=7.3, height = 4.3)

svz3 <- FindNeighbors(svz3, dims = 1:20)
svz3 <- FindClusters(svz3, resolution = 0.15)
svz3 <- RunUMAP(svz3, dims = 1:20)
DimPlot(svz3)
ggsave("plots1/umap.initial.fil.lane3.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz3), colnames(lmo.data3)) #1462
lmo.data3 <- lmo.data3[ , barcode.intersect] 
svz3[["LMO"]] <- CreateAssayObject(counts = lmo.data3)

#MULTIseq
svz3 <- NormalizeData(object = svz3, assay = "LMO", normalization.method = "CLR", margin=2)

multi3 <- MULTIseqDemux(svz3, assay = "LMO", 
                        autoThresh = T, maxiter = 2, 
                        qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi3$MULTI_ID)
# D-L0112-28.1-M-lane3-Dox-ACCAATGC  Doublet 
# 48                                13 
# F-L0116-28.1-F-lane3-Dox-CGAACAAG  H-L125-27.4-M-lane3-Dox-GAAGCTTG 
# 107                               143 
# Negative 
# 1151 

# Semi-supervised negative barcode reclassification
df <- t(as.matrix(multi3@assays$LMO@data))
bar.table.full <- df[,1:3]
df2 <- as.data.frame(multi3$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi3$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 6) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]


multi3$MULTI_classification_rescued <- final.calls.rescued
table(multi3$MULTI_classification_rescued)
# D-L0112-28.1-M-lane3-Dox-ACCAATGC                           Doublet 
# 285                                13 
# F-L0116-28.1-F-lane3-Dox-CGAACAAG  H-L125-27.4-M-lane3-Dox-GAAGCTTG 
# 275                               210 
# Negative 
# 679 
table(multi3$MULTI_classification)

Idents(multi3) <- "MULTI_classification_rescued"

multi3$MULTI_classification.global <- plyr::mapvalues(multi3$MULTI_classification_rescued,
                                                      from = c("D-L0112-28.1-M-lane3-Dox-ACCAATGC",
                                                               "F-L0116-28.1-F-lane3-Dox-CGAACAAG",
                                                               "H-L125-27.4-M-lane3-Dox-GAAGCTTG",
                                                               "Negative","Doublet"),
                                                      to = c("Singlet", "Singlet", "Singlet",
                                                             "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi3$MULTI_classification.global <- factor(multi3$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi3) <- "MULTI_classification.global"
VlnPlot(multi3, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

VlnPlot(multi3, features = "nFeature_RNA", pt.size = 0.1, log = F)

VlnPlot(multi3, features = "nFeature_RNA", pt.size = 0.1, log = T)


#===================================
# Intermediate Saving/Loading Point
saveRDS(svz3, "data/svz3.05192022.rds")

saveRDS(multi3, "data/svz3.MULTI.CLR.2.rescued.05192022.rds")


#=======================================================================
#=======================================================================
# 10x Lane 4 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data4 <- Read10X("../1.CellRanger/output4/filtered_feature_bc_matrix") #31059 10165, genes x cells

lmo.data4 <- Read10X("../2.LMO/output4/umi_count/", gene.column=1) #4 10165, samples and unmapped x cells
lmo.data4 <- lmo.data4[1:3,] #3 10165, samples x cells

colnames(svz.data4) <- sub("-1","",colnames(svz.data4))
barcode.intersect <- intersect(colnames(svz.data4), colnames(lmo.data4)) #10165
lmo.data4 <- lmo.data4[ , barcode.intersect] 
svz.data4 <- svz.data4[ , barcode.intersect] 

svz4 <- CreateSeuratObject(counts = svz.data4, project = "Lane4") #31059 features across 10165 samples within 1 assay
svz4 <- NormalizeData(svz4)
svz4 <- FindVariableFeatures(svz4)
svz4 <- ScaleData(svz4)


#======================================================================================
# QC Metrics
#======================================================================================

svz4[["percent.mt"]] <- PercentageFeatureSet(svz4, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(svz4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.lane4.pdf", width = 7, height = 5)

# cut off of less than 10% mito; more than 500 features.
svz4 <- subset(svz4, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 9582 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane4.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select PCs for tSNE visualization and graph-based clustering
svz4 <- RunPCA(svz4, verbose = FALSE)
ElbowPlot(svz4, ndims = 50)
ggsave("plots1/elbow.lane4.pdf", width=7.3, height = 4.3)

svz4 <- FindNeighbors(svz4, dims = 1:19)
svz4 <- FindClusters(svz4, resolution = 0.15)
svz4 <- RunUMAP(svz4, dims = 1:19)
DimPlot(svz4)
ggsave("plots1/umap.initial.fil.lane4.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz4), colnames(lmo.data4)) #9582
lmo.data4 <- lmo.data4[ , barcode.intersect] 
svz4[["LMO"]] <- CreateAssayObject(counts = lmo.data4)

#MULTIseq----
svz4 <- NormalizeData(object = svz4, assay = "LMO", normalization.method = "CLR", margin=2)

multi4 <- MULTIseqDemux(svz4, assay = "LMO", 
                        autoThresh = T, maxiter = 5, 
                        qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi4$MULTI_ID)
# Doublet                          I-1114L-4.0-F-lane4-untr-AAGTACGC 
# 90                               608 
# J-1115L-4.0-F-lane4-untr-ATTCGCAC K-1116L-4.0-F-lane4-untr-GAGTCGAT 
# 580                               466 
# Negative 
# 7838 

# Semi-supervised negative barcode reclassification
df <- t(as.matrix(multi4@assays$LMO@data))
bar.table.full <- df[,1:3]
df2 <- as.data.frame(multi4$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi4$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 4) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

multi4$MULTI_classification_rescued <- final.calls.rescued
table(multi4$MULTI_classification_rescued)
# Doublet                         I-1114L-4.0-F-lane4-untr-AAGTACGC 
# 90                              2441 
# J-1115L-4.0-F-lane4-untr-ATTCGCAC K-1116L-4.0-F-lane4-untr-GAGTCGAT 
# 2412                              2212 
# Negative 
# 2427 

Idents(multi4) <- "MULTI_classification_rescued"

multi4$MULTI_classification.global <- plyr::mapvalues(multi4$MULTI_classification_rescued,
                                                      from = c("I-1114L-4.0-F-lane4-untr-AAGTACGC",
                                                               "J-1115L-4.0-F-lane4-untr-ATTCGCAC",
                                                               "K-1116L-4.0-F-lane4-untr-GAGTCGAT",
                                                               "Negative","Doublet"),
                                                      to = c("Singlet", "Singlet", "Singlet",
                                                             "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi4$MULTI_classification.global <- factor(multi4$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi4) <- "MULTI_classification.global"
VlnPlot(multi4, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

VlnPlot(multi4, features = "nFeature_RNA", pt.size = 0.1, log = F)

VlnPlot(multi4, features = "nFeature_RNA", pt.size = 0.1, log = T)

#===================================
# Intermediate Saving/Loading Point
saveRDS(svz4, "data/svz4.05192022.rds")

saveRDS(multi4, "data/svz4.MULTI.CLR.2.rescued.05192022.rds")


#=======================================================================
#=======================================================================
# 10x Lane 5 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data5 <- Read10X("../1.CellRanger/output5/filtered_feature_bc_matrix") #31059  2715, genes x cells

lmo.data5 <- Read10X("../2.LMO/output5/umi_count/", gene.column=1) #4 2715, samples and unmapped x cells
lmo.data5 <- lmo.data5[1:3,] #3 2715, samples x cells

colnames(svz.data5) <- sub("-1","",colnames(svz.data5))
barcode.intersect <- intersect(colnames(svz.data5), colnames(lmo.data5)) 
lmo.data5 <- lmo.data5[ , barcode.intersect] 
svz.data5 <- svz.data5[ , barcode.intersect] 

svz5 <- CreateSeuratObject(counts = svz.data5, project = "Lane5") #31059 features across 2715 samples within 1 assay
svz5 <- NormalizeData(svz5)
svz5 <- FindVariableFeatures(svz5)
svz5 <- ScaleData(svz5)


#======================================================================================
# QC Metrics
#======================================================================================

svz5[["percent.mt"]] <- PercentageFeatureSet(svz5, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(svz5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.lane5.pdf", width = 7, height = 5)

# cut off of less than 10% mito; more than 500 features.
svz5 <- subset(svz5, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 2240 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane5.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 18 PCs for tSNE visualization and graph-based clustering
svz5 <- RunPCA(svz5, verbose = FALSE)
ElbowPlot(svz5, ndims = 50)
ggsave("plots1/elbow.lane5.pdf", width=7.3, height = 4.3)

svz5 <- FindNeighbors(svz5, dims = 1:18)
svz5 <- FindClusters(svz5, resolution = 0.15)
svz5 <- RunUMAP(svz5, dims = 1:18)
DimPlot(svz5)
ggsave("plots1/umap.initial.fil.lane5.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz5), colnames(lmo.data5)) #2240
lmo.data5 <- lmo.data5[ , barcode.intersect] 
svz5[["LMO"]] <- CreateAssayObject(counts = lmo.data5)

#MULTIseq
svz5 <- NormalizeData(object = svz5, assay = "LMO", normalization.method = "CLR", margin=2)

multi5 <- MULTIseqDemux(svz5, assay = "LMO", 
                        autoThresh = T, maxiter = 5, 
                        qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi5$MULTI_ID)
# Doublet  M-L126-27.4-M-lane5-untr-CAGTTAGG 
# 57                                209 
# Negative  O-L132-27.4-F-lane5-untr-AAGCAGTC 
# 1600                                323 
# Q-L0108-28.8-M-lane5-untr-ACTCGAAG 
# 51 

Idents(multi5) <- "MULTI_ID"
RidgePlot(multi5, assay = "LMO", features = rownames(multi5[["LMO"]])[1:3], ncol = 1)
ggsave("plots1/MULTI.CLR.2/ridge.multi.CLR.lane5.pdf", height = 7, width = 10)

DoHeatmap(multi5, assay = "LMO", slot = "data", group.by="MULTI_ID", 
          features=rownames(multi5[["LMO"]])[1:3]
          , draw.lines = F, size=3, angle=20) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots1/MULTI.CLR.2/heatmap.multi.lane5.png", width = 10.9, height=4)

# Semi-supervised negative barcode reclassification
df <- t(as.matrix(multi5@assays$LMO@data))
bar.table.full <- df[,1:3]
df2 <- as.data.frame(multi5$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi5$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 5) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

multi5$MULTI_classification_rescued <- final.calls.rescued
table(multi5$MULTI_classification_rescued)
# Doublet  M-L126-27.4-M-lane5-untr-CAGTTAGG 
# 57                                480 
# Negative  O-L132-27.4-F-lane5-untr-AAGCAGTC 
# 402                                953 
# Q-L0108-28.8-M-lane5-untr-ACTCGAAG 
# 348 

Idents(multi5) <- "MULTI_classification_rescued"

multi5$MULTI_classification.global <- plyr::mapvalues(multi5$MULTI_classification_rescued,
                                                      from = c("M-L126-27.4-M-lane5-untr-CAGTTAGG",
                                                               "O-L132-27.4-F-lane5-untr-AAGCAGTC",
                                                               "Q-L0108-28.8-M-lane5-untr-ACTCGAAG",
                                                               "Negative","Doublet"),
                                                      to = c("Singlet", "Singlet", "Singlet",
                                                             "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi5$MULTI_classification.global <- factor(multi5$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi5) <- "MULTI_classification.global"
VlnPlot(multi5, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

VlnPlot(multi5, features = "nFeature_RNA", pt.size = 0.1, log = F)

VlnPlot(multi5, features = "nFeature_RNA", pt.size = 0.1, log = T)



#===================================
# Intermediate Saving/Loading Point
saveRDS(svz5, "data/svz5.05192022.rds")

saveRDS(multi5, "data/svz5.MULTI.CLR.2.rescued.05192022.rds")


#=======================================================================
#=======================================================================
# 10x Lane 6 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data6 <- Read10X("../1.CellRanger/output6/filtered_feature_bc_matrix") #31059  3159, genes x cells

lmo.data6 <- Read10X("../2.LMO/output6/umi_count/", gene.column=1) #4 3159, samples and unmapped x cells
lmo.data6 <- lmo.data6[1:3,] #3 3159, samples x cells

colnames(svz.data6) <- sub("-1","",colnames(svz.data6))
barcode.intersect <- intersect(colnames(svz.data6), colnames(lmo.data6)) 
lmo.data6 <- lmo.data6[ , barcode.intersect] 
svz.data6 <- svz.data6[ , barcode.intersect] 

svz6 <- CreateSeuratObject(counts = svz.data6, project = "Lane6") #31059 features across 3159 samples within 1 assay
svz6 <- NormalizeData(svz6)
svz6 <- FindVariableFeatures(svz6)
svz6 <- ScaleData(svz6)


#======================================================================================
# QC Metrics
#======================================================================================

svz6[["percent.mt"]] <- PercentageFeatureSet(svz6, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(svz6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.lane6.pdf", width = 7, height = 5)

# cut off of less than 10% mito; more than 500 features.
svz6 <- subset(svz6, subset = percent.mt < 10 & nFeature_RNA > 500) 
#31059 features across 2992 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots1/violin.mt.feature.counts.fil.lane6.pdf", width = 7, height = 5)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 19 PCs for tSNE visualization and graph-based clustering
svz6 <- RunPCA(svz6, verbose = FALSE)
ElbowPlot(svz6, ndims = 50)
ggsave("plots1/elbow.lane6.pdf", width=7.3, height = 4.3)

svz6 <- FindNeighbors(svz6, dims = 1:19)
svz6 <- FindClusters(svz6, resolution = 0.15)
svz6 <- RunUMAP(svz6, dims = 1:19)
DimPlot(svz6)
ggsave("plots1/umap.initial.fil.lane6.pdf", width = 5, height = 5)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz6), colnames(lmo.data6)) #2992
lmo.data6 <- lmo.data6[ , barcode.intersect] 
svz6[["LMO"]] <- CreateAssayObject(counts = lmo.data6)

#MULTIseq
svz6 <- NormalizeData(object = svz6, assay = "LMO", normalization.method = "CLR", margin=2)

multi6 <- MULTIseqDemux(svz6, assay = "LMO", 
                        autoThresh = T, maxiter = 50, 
                        qrange = seq(from = 0.1, to = 0.9, by = 0.05), verbose = TRUE)

table(multi6$MULTI_ID)
# Doublet L-L131-27.4-F-lane6-Dox-AAGGCTAG 
# 3                               40 
# N-L128-27.4-F-lane6-Dox-AACCGAAC                         Negative 
# 7                             2933 
# P-L133-26.3-M-lane6-Dox-GAATCAGG 
# 9 

# Semi-supervised negative barcode reclassification
df <- t(as.matrix(multi6@assays$LMO@data))
bar.table.full <- df[,1:3]
df2 <- as.data.frame(multi6$MULTI_ID)
df2 <- tibble::rownames_to_column(df2)
final.calls <- df2$`multi6$MULTI_ID`
names(final.calls) <- df2$rowname

reclass.cells <- findReclassCells(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

str(reclass.cells)
str(reclass.res)

## Visualize Results
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
  geom_point() + 
  xlim(c(nrow(reclass.res)-1,1)) + 
  # ylim(c(0,1.05)) +
  geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
  geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)

## Finalize negative cell rescue results
final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 5) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]

multi6$MULTI_classification_rescued <- final.calls.rescued
table(multi6$MULTI_classification_rescued)
# Doublet L-L131-27.4-F-lane6-Dox-AAGGCTAG 
# 3                              666 
# N-L128-27.4-F-lane6-Dox-AACCGAAC                         Negative 
# 479                              495 
# P-L133-26.3-M-lane6-Dox-GAATCAGG 
# 1349 

Idents(multi6) <- "MULTI_classification_rescued"

multi6$MULTI_classification.global <- plyr::mapvalues(multi6$MULTI_classification_rescued,
                                                      from = c("L-L131-27.4-F-lane6-Dox-AAGGCTAG",
                                                               "N-L128-27.4-F-lane6-Dox-AACCGAAC",
                                                               "P-L133-26.3-M-lane6-Dox-GAATCAGG",
                                                               "Negative","Doublet"),
                                                      to = c("Singlet", "Singlet", "Singlet",
                                                             "Negative","Doublet"))

# Violin plots of expression by single/doublet/negative catagories
multi6$MULTI_classification.global <- factor(multi6$MULTI_classification.global,order=T,levels=c("Singlet","Doublet","Negative"))
Idents(multi6) <- "MULTI_classification.global"
VlnPlot(multi6, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

VlnPlot(multi6, features = "nFeature_RNA", pt.size = 0.1, log = F)

VlnPlot(multi6, features = "nFeature_RNA", pt.size = 0.1, log = T)


#===================================
# Intermediate Saving/Loading Point
saveRDS(svz6, "data/svz6.05192022.rds")

saveRDS(multi6, "data/svz6.MULTI.CLR.2.rescued.05192022.rds")


#======================================================================================
# Load Seurat objects
#======================================================================================

svz1 <- readRDS("data/svz1.MULTI.CLR.2.rescued.05192022.rds")
svz2 <- readRDS("data/svz2.MULTI.CLR.2.rescued.05192022.rds")
svz3 <- readRDS("data/svz3.MULTI.CLR.2.rescued.05192022.rds")
svz4 <- readRDS("data/svz4.MULTI.CLR.2.rescued.05192022.rds")
svz5 <- readRDS("data/svz5.MULTI.CLR.2.rescued.05192022.rds")
svz6 <- readRDS("data/svz6.MULTI.CLR.2.rescued.05192022.rds")



#======================================================================================
# Merge Seurat objects
#======================================================================================

svz <- merge(svz1, c(svz2,svz3,svz4,svz5,svz6), add.cell.ids = c("l1","l2","l3","l4","l5","l6"), project = "SVZ.OSKM", merge.data = TRUE)
#31076 features across 22080 samples within 2 assays 

# standard log-normalization
svz <- NormalizeData(svz)
# choose ~1k variable features
svz <- FindVariableFeatures(svz)
# standard scaling (no regression)
svz <- ScaleData(svz)

# Run PCA
svz <- RunPCA(svz, verbose = FALSE)
ElbowPlot(svz, ndims = 50)
ggsave("plots1/elbow.combined.pdf", width=7.3, height = 4.3)

svz <- FindNeighbors(svz, dims = 1:23)
svz <- FindClusters(svz, resolution = 0.15)
svz <- RunUMAP(svz, dims = 1:23)
DimPlot(svz)
ggsave("plots1/umap.initial.fil.combined.pdf", width = 7, height = 7)


#======================================================================================

# UMAP based on transcriptome reduction, colored by global multiseq classification
Idents(svz) <- "MULTI_classification.global"
DimPlot(svz, cols = colors3, reduction = "umap")
ggsave("plots1/umap.doub.sing.neg.combined.pdf", width = 6.5, height = 6)

Idents(svz) <- "MULTI_classification.global"
DimPlot(svz, cols = colors3, reduction = "umap", split.by = "MULTI_classification.global")
ggsave("plots1/umap.split.combined.pdf", width = 7, height = 3)

saveRDS(svz, paste0("data/svz.combined.", Sys.Date(), ".rds"))


#======================================================================================
# Subset to singlets
#======================================================================================

# Set colors
color_pal.1 <- tableau_color_pal(palette = "Tableau 10")(1)
color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Remove doublets and Negatives
svz_sing <- subset(svz, subset = MULTI_classification.global == "Singlet") 
#31076 features across 16291 samples within 2 assays

# UMAP based on transcriptome reduction, colored by doublet/singlet status
DimPlot(svz_sing, cols = colors3, reduction = "umap", group.by="MULTI_classification.global")
ggsave("plots1/umap.sing.combined.pdf", width = 5.9, height = 4.6)

# UMAP based on transcriptome reduction, colored by sample barcode
DimPlot(svz_sing, reduction = "umap", group.by="MULTI_classification_rescued")
ggsave("plots1/umap.sing.combined.sample.pdf", width = 9.2, height = 4.6)

# UMAP based on transcriptome reduction, split by lane
DimPlot(svz_sing, reduction = "umap", split.by="orig.ident", group.by="MULTI_classification.global", cols = colors3)
ggsave("plots1/umap.sing.combined.lanesplit.pdf", width = 9.2, height = 3)


table(svz_sing$MULTI_classification_rescued)
# A-1108L-3.9-M-lane1-untr-TGTGATGG  B-1110L-3.9-M-lane1-untr-TCAATGGC 
# 1422                               1428 
# C-L0111-28.1-M-lane2-untr-CTCTAGAC  D-L0112-28.1-M-lane3-Dox-ACCAATGC 
# 712                                285 
# E-L0114-28.1-M-lane2-untr-AGTTGCGT  F-L0116-28.1-F-lane3-Dox-CGAACAAG 
# 209                                275 
# G-L0117-28.1-F-lane2-untr-GTACCTGT   H-L125-27.4-M-lane3-Dox-GAAGCTTG 
# 410                                210 
# I-1114L-4.0-F-lane4-untr-AAGTACGC  J-1115L-4.0-F-lane4-untr-ATTCGCAC 
# 2441                               2412 
# K-1116L-4.0-F-lane4-untr-GAGTCGAT   L-L131-27.4-F-lane6-Dox-AAGGCTAG 
# 2212                                666 
# M-L126-27.4-M-lane5-untr-CAGTTAGG   N-L128-27.4-F-lane6-Dox-AACCGAAC 
# 480                                479 
# O-L132-27.4-F-lane5-untr-AAGCAGTC   P-L133-26.3-M-lane6-Dox-GAATCAGG 
# 953                               1349 
# Q-L0108-28.8-M-lane5-untr-ACTCGAAG 
# 348 



#======================================================================================
# Rerun processing on subset
#======================================================================================

# Visualize QC metrics as a violin plot
VlnPlot(svz_sing, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=.1)
ggsave("plots1/violin.mt.feature.counts.sing.sample.pdf", width = 20, height=10)

# Normalize, Standardize. SCTransform will help with different sequencing depths. 
svz_sing <- SCTransform(svz_sing)

# Reduce
svz_sing <- RunPCA(svz_sing, verbose = FALSE)
ElbowPlot(svz_sing, ndims = 50)
ggsave("plots1/elbow.sing.combined.sct.pdf", width = 5.3, height = 2.8)

svz_sing <- FindNeighbors(svz_sing, dims = 1:14)
svz_sing <- FindClusters(svz_sing, resolution = 0.5)
svz_sing <- RunUMAP(svz_sing, dims = 1:14)

DimPlot(svz_sing, pt.size = .6)
ggsave("plots1/umap.sing.combined.sct.pdf", width = 5.9, height = 4.6)

DimPlot(svz_sing, pt.size = .6, label=T)
ggsave("plots1/umap.sing.combined.sct.numbered.pdf", width = 5.9, height = 4.6)

DimPlot(svz_sing, group.by = "MULTI_classification_rescued", pt.size = .6)
ggsave("plots1/umap.sing.combined.sct.sample.pdf", width = 10, height = 6)

DimPlot(svz_sing, group.by = "orig.ident", pt.size = .6)
ggsave("plots1/umap.sing.combined.sct.lane.pdf", width = 5.9, height = 4.6)

DimPlot(svz_sing, split.by = "orig.ident", pt.size = .6)
ggsave("plots1/umap.sing.combined.sct.lanesplit.pdf", width = 12, height = 4.6)


#-----
svz_sing.markers <- FindAllMarkers(object=svz_sing)

saveRDS(svz_sing.markers, paste0("data/svz.combined.sing.MULTI.CLR.2.rescued.markers.", Sys.Date(), ".rds"))
write.csv(svz_sing.markers, paste0('data/svz.combined.sing.MULTI.CLR.2.rescued.markers.', Sys.Date(), '.csv'))

# Pull out top 10 for heatmap 
top10 <- svz_sing.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(svz_sing, features = top10$gene, lines.width = 10) + NoLegend() + theme(axis.text.y = element_text(size = 10))
ggsave("plots1/heatmap.markers.top10.label.png", height = 20, width = 20)
write.csv(top10, paste0('data/top10_svz.combined.sing.MULTI.CLR.2.rescued.markers.', Sys.Date(), '.csv'))


#======================================================================================
# Save Seurat object
#======================================================================================
saveRDS(svz_sing, paste0("data/svz.combined.sing.MULTI.CLR.2.rescued.", Sys.Date(), ".rds"))
#svz_sing <- readRDS("data/svz.combined.sing.2021-08-30.rds")


#===============
sessionInfo()
#  R version 4.0.2 (2020-06-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# [1] cowplot_1.1.0               viridis_0.5.1               viridisLite_0.3.0          
# [4] ggthemes_4.2.0              scales_1.1.1                sctransform_0.3.1          
# [7] Matrix_1.2-18               forcats_0.5.0               stringr_1.4.0              
# [10] purrr_0.3.4                 readr_1.3.1                 tidyr_1.1.2                
# [13] tibble_3.0.4                ggplot2_3.3.2               tidyverse_1.3.0            
# [16] MAST_1.14.0                 SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.2
# [19] DelayedArray_0.14.1         matrixStats_0.57.0          Biobase_2.48.0             
# [22] GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
# [25] S4Vectors_0.26.1            BiocGenerics_0.34.0         dplyr_1.0.2                
# [28] Seurat_3.2.2                 deMULTIplex_1.0.2
