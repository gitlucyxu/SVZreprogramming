# Whole-body (iOSKM) partial reprogramming, cohort 1
# LMO + 6 SVZ Seurat Processing

# VERSION 3.1.2 Seurat

rm(list = ls())
library(Seurat)
library(dplyr)
library(MAST)
library(tidyverse)
library(Matrix)
library(sctransform)
library(scales)
library(ggthemes)
library(viridis)
library(cowplot)


# All subsequent paths will be relative.
setwd("~/Dropbox/OSKM_10x/3.Seurat")


#======================================================================================
# Create Seurat object and basic processing
#======================================================================================
# Using Cellranger run that mapped to custom genome (mm10, FUW-M2rtTA full-length, and each 2A sequence)
# Using LMO counts from "Undetermined" files in GEX library
#======

# Read in Cellranger count matrix and create basic seurat object
svz.data <- Read10X("../1.CellRanger/output_04012020") #31054 x 10842, genes x cells

# Read in LMO count matrix, remove unmapped row
lmo.data2 <- Read10X("../2.LMO/output2/umi_count/", gene.column=1) #7 x 10835, samples and unmapped x cells
lmo.data2 <- lmo.data2[1:6,] #6 x 10830, samples x cells

# Remove cells without barcodes and barcodes without cells
barcode.intersect <- intersect(colnames(svz.data), colnames(lmo.data2))
lmo.data2 <- lmo.data2[ , barcode.intersect] #6 x 10830
svz.data <- svz.data[ , barcode.intersect] #31054 x 10830

# Make Seurat object
svz <- CreateSeuratObject(counts = svz.data) #31054 features across 10830 samples within 1 assay
svz <- NormalizeData(svz) # standard log-normalization
svz <- FindVariableFeatures(svz) # choose ~1k variable features
svz <- ScaleData(svz) # standard scaling (no regression)


#======================================================================================
# QC Metrics and Filtering
#======================================================================================

svz[["percent.mt"]] <- PercentageFeatureSet(svz, pattern = "^mt-")

# Cut off of 1000 UMI; less than 10% mito; more than 500 features.
svz <- subset(svz, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 
#31057 features across 9284 samples within 1 assay


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 17 PCs for tSNE visualization and graph-based clustering
svz <- RunPCA(svz, verbose = FALSE)
ElbowPlot(svz, ndims = 50)

svz <- FindNeighbors(svz, dims = 1:17)
svz <- FindClusters(svz, resolution = 0.15)
svz <- RunUMAP(svz, dims = 1:17)
DimPlot(svz)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz), colnames(lmo.data2)) #9284
lmo.data2 <- lmo.data2[ , barcode.intersect] #6 x 9284
svz[["LMO2"]] <- CreateAssayObject(counts = lmo.data2)
#An object of class Seurat 
#31063 features across 9284 samples within 2 assays 
#Active assay: RNA (31057 features)
#1 other assay present: LMO2
#2 dimensional reductions calculated: pca, umap


svz <- NormalizeData(object = svz, assay = "LMO2", normalization.method = "CLR")
svz <- HTODemux(svz, assay = "LMO2", positive.quantile = 0.8)
# Doublet Negative  Singlet 
# 4979     2477     1828  

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz, assay = "LMO2", ncells = 10000) + scale_fill_viridis(option="viridis", direction=1)

# Some droplets have very high expression of all LMOs- subset based on nCount_LMO2 < 2000 to try to get rid of these before determining threshold for positive LMO reads
svz.subset.lmocount2000 <- subset(svz, subset = nCount_LMO2 < 2000)
# 31063 features across 5323 samples within 2 assays 

# Re-normalize RNA
svz.subset.lmocount2000 <- NormalizeData(svz.subset.lmocount2000)
svz.subset.lmocount2000 <- FindVariableFeatures(svz.subset.lmocount2000)
svz.subset.lmocount2000 <- ScaleData(svz.subset.lmocount2000)

# Run PCA, select 16 PCs for tSNE visualization and graph-based clustering
svz.subset.lmocount2000 <- RunPCA(svz.subset.lmocount2000, verbose = FALSE)
ElbowPlot(svz.subset.lmocount2000, ndims = 50)

svz.subset.lmocount2000 <- FindNeighbors(svz.subset.lmocount2000, dims = 1:16)
svz.subset.lmocount2000 <- FindClusters(svz.subset.lmocount2000, resolution = 0.15)
svz.subset.lmocount2000 <- RunUMAP(svz.subset.lmocount2000, dims = 1:16)
DimPlot(svz.subset.lmocount2000)

#Re-normalize LMOs
svz.subset.lmocount2000 <- NormalizeData(object = svz.subset.lmocount2000, assay = "LMO2", normalization.method = "CLR")
svz.subset.lmocount2000 <- HTODemux(svz.subset.lmocount2000, assay = "LMO2", positive.quantile = 0.95)
# Doublet Negative  Singlet 
# 1402     1671     2250 

# Heatmap of LMO expression
Idents(svz.subset.lmocount2000) <- "LMO2_maxID"
HTOHeatmap(svz.subset.lmocount2000, assay = "LMO2", ncells = 10000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots_04012020/subsetlmocount2000.heatmap.lmo.all.png", width = 10.9, height=2.88)


# Violin plots of expression by single/doublet/negative categories
Idents(svz.subset.lmocount2000) <- "LMO2_classification.global"
VlnPlot(svz.subset.lmocount2000, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


#======================================================================================
# Reductions based on LMO counts

# Set three colors
color_pal.3 <- tableau_color_pal(palette = "Tableau 10")(3)[c(2,3,1)]

# Calculate a distance matrix using HTO
lmo2.dist.mtx <- as.matrix(dist(t(GetAssayData(object = svz.subset.lmocount2000, assay = "LMO2"))))

# Calculate tSNE embeddings with distance matrix
svz.subset.lmocount2000 <- RunTSNE(svz.subset.lmocount2000, distance.matrix = lmo2.dist.mtx, perplexity = 50)

# tSNE based on LMO reduction, colored by global lmo classifiction
DimPlot(svz.subset.lmocount2000, cols = color_pal.3, reduction = "tsne")

# UMAP based on transcriptome reduction, colored by global lmo classification
DimPlot(svz.subset.lmocount2000, cols = color_pal.3, reduction = "umap")


#======================================================================================
# Subset to singlets
#======================================================================================

# Set colors
color_pal.1 <- tableau_color_pal(palette = "Tableau 10")(1)
color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Remove doublets and Negatives
svz_sing <- subset(svz.subset.lmocount2000, subset = LMO2_classification.global == "Singlet") 
#31063 features across 2250 samples within 2 assays 

# Calculate a distance matrix using HTO
lmo.dist.mtx.sing <- as.matrix(dist(t(GetAssayData(object = svz_sing, assay = "LMO2"))))

# Calculate tSNE embeddings with a distance matrix
svz_sing <- RunTSNE(svz_sing, distance.matrix = lmo.dist.mtx.sing, perplexity = 50)

# tSNE based on LMO reduction, colored by doublet/singlet status
DimPlot(svz_sing, cols = color_pal.1, reduction = "tsne")

# tSNE based on LMO reduction, colored by sample barcode
DimPlot(svz_sing, cols = color_pal.6, reduction = "tsne", group.by="hash.ID")

# UMAP based on transcriptome reduction, colored by doublet/singlet status
DimPlot(svz_sing, cols = color_pal.1, reduction = "umap")

# UMAP based on transcriptome reduction, colored by sample barcode
DimPlot(svz_sing, cols = color_pal.6, reduction = "umap", group.by="hash.ID")


table(svz_sing$LMO2_maxID)
# Sample1-961L-DOB20180508-2Dox0-CTCTAGAC  Sample2-955L-DOB20180501-untr-ACCAATGC Sample3-928L-DOB20180327-2Dox5-AGTTGCGT 
# 212                                     347                                     515 
# Sample4-960L-DOB20180508-2Dox0-CGAACAAG Sample5-929L-DOB20180327-2Dox5-GTACCTGT  Sample6-953L-DOB20180501-untr-GAAGCTTG 
# 432                                     541                                     203 

table(svz_sing$Age_Treatment)
#Total untr     Total 2Dox0       Total 2Dox5
#550            644               1056



#======================================================================================
# Rerun processing on fully filtered subset
#======================================================================================

# Normalize, Standardize
svz_sing <- SCTransform(svz_sing)

# Reduce
svz_sing <- RunPCA(svz_sing, verbose = FALSE)
ElbowPlot(svz_sing, ndims = 50)

svz_sing <- FindNeighbors(svz_sing, dims = 1:18)
svz_sing <- FindClusters(svz_sing, resolution = 0.5)
svz_sing <- RunUMAP(svz_sing, dims = 1:18)

DimPlot(svz_sing, pt.size = .6)
DimPlot(svz_sing, group.by = "hash.ID", pt.size = .6, cols = color_pal.6)


#======================================================================================
# Adjust settings of UMAP to separate small clusters

svz_sing <- FindNeighbors(svz_sing, dims = 1:18)
svz_sing_0.8 <- FindClusters(svz_sing, resolution = 0.8)
svz_sing_0.8 <- RunUMAP(svz_sing_0.8, dims = 1:18, min.dist = .8, spread = .5, seed.use = 55)

DimPlot(svz_sing_0.8, pt.size = .6, label = T)

DimPlot(svz_sing_0.8, group.by = "hash.ID", pt.size = .6, cols = color_pal.6)

#======================================================================================
# Save Seurat object
#======================================================================================
saveRDS(svz_sing_0.8, paste0("data_04012020/svz_sing.0.8_", Sys.Date(), ".rds"))
svz_sing_0.8 <- readRDS("data_04012020/svz_sing.0.8_2020-04-12.rds")

#======================================================================================
# Find markers
#======================================================================================

# Find marker genes for clusters and save
svz_sing.0.8.markers <- FindAllMarkers(object=svz_sing_0.8)
write.csv(svz_sing.0.8.markers, paste0('data_04012020/svz_sing.0.8.markers.', Sys.Date(), '.csv'))
saveRDS(svz_sing.0.8.markers, paste0("data_04012020/svz_sing.0.8.markers_", Sys.Date(), ".rds"))

# Pull out top 10 for heatmap 
top10_0.8 <- svz_sing.0.8.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz_sing_0.8, features = top10_0.8$gene, lines.width = 10, label = F) + NoLegend() + theme(axis.text.y = element_text(size = 10))
ggsave("plots_04012020/heatmap.0.8.markers.png", height = 20, width = 20)
write.csv(top10_0.8, paste0('data_04012020/top10_svz_sing.0.8.markers.', Sys.Date(), '.csv'))

# Pull out top 5 for heatmap 
top5 <- svz_sing.0.8.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(svz_sing_0.8, features = top5$gene, lines.width = 10, label=F) + NoLegend() + theme(axis.text.y = element_text(size = 14))
ggsave("plots_04012020/heatmap.0.8.markers.top5.png", height = 20, width = 20)



#===== Cell type markers =====#

# Microglia, olig, and endothelial marker genes
genes <- c("C1qc", "Junb", "Cd52", "Opalin", "Tspan2", "Klk6", "Cxcl12", "Ly6c1", "Itm2a")
FeaturePlot(svz, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.all.micro-olig-endo-markers.pdf", width = 12, height = 10)

genes <- c("C1qc", "Opalin", "Ly6c1")
FeaturePlot(svz.subset.lmocount2000, features = genes, split.by = "LMO2_classification.global", order = TRUE, pt.size = 0.8)
ggsave("plots2_04012020/umap.subset-doub-neg-sing.C1qc-Opalin-Ly6c1.pdf", width = 12, height = 12)


# NSC lineage 
genes <- c("Clu", "Aldoc", "Meg3", "Tubb2b", "Pclaf", "Cenpf")
FeaturePlot(svz, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.all.NSClineage-markers.pdf", width = 8, height = 10)
FeaturePlot(svz.subset.lmocount2000, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.subset.NSClineage-markers.pdf", width = 8, height = 10)


# Small clusters: mural, neurons?, OPCs?
genes <- c("Myl9", "Acta2", "Vtn", "Rarres2", "Ttr", "1500015O10Rik", "Gpr17", "Tnr", "Ptprz1")
FeaturePlot(svz, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.all.mural-neuron-OPC-markers.pdf", width = 12, height = 10)

genes <- c("Myl9", "1500015O10Rik", "Tnr")
FeaturePlot(svz.subset.lmocount2000, features = genes, split.by = "LMO2_classification.global", order = TRUE, pt.size = 0.8)
ggsave("plots2_04012020/umap.subset-doub-neg-sing.Myl9-1500015O10Rik-Tnr.pdf", width = 12, height = 12)


# T cells
genes <- c("Cd3g", "Hcst", "Ccl5", "Nkg7", "AW112010")
FeaturePlot(svz, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.all.Tcell-markers.pdf", width = 8, height = 10)

# T cell genes: Cd3, Cd8, Cd4, interferon, activation, and tissue retention genes
cd.genes <- grep(pattern = "^Cd3", x = rownames(x = svz$SCT@data), value = TRUE)
genes <- c("Cd3g", "Cd3d", "Cd3e", "Cd8a", "Cd4", "Ifng", "Cd69", "Xcl1", "Itga4") #"Cd8b"
FeaturePlot(svz, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.all.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)
FeaturePlot(svz.subset.lmocount2000, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.subset.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)
FeaturePlot(svz_sing_0.8, features = genes, order = TRUE)
ggsave("plots2_04012020/umap.sing.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)





#----
sessionInfo()

#  R version 3.5.0 (2018-04-23)
# [1] viridis_0.5.1               viridisLite_0.3.0           ggthemes_4.2.0              scales_1.1.0               
# [5] sctransform_0.2.1           Matrix_1.2-18               MAST_1.8.2                  SingleCellExperiment_1.4.1 
# [9] SummarizedExperiment_1.10.1 DelayedArray_0.6.1          BiocParallel_1.14.1         matrixStats_0.53.1         
# [13] Biobase_2.40.0              GenomicRanges_1.32.3        GenomeInfoDb_1.16.0         IRanges_2.14.10            
# [17] S4Vectors_0.18.3            BiocGenerics_0.26.0         gridExtra_2.3               ggloop_0.1.0               
# [21] cowplot_0.9.2               forcats_0.4.0               stringr_1.4.0               dplyr_0.8.3                
# [25] purrr_0.3.3                 readr_1.3.1                 tidyr_1.0.0                 tibble_2.1.3               
# [29] ggplot2_3.3.0               tidyverse_1.3.0             Seurat_3.1.2 
