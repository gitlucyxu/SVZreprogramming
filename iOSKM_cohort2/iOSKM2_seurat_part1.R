# Whole-body (iOSKM) partial reprogramming, cohort 2
# LMO + 12 SVZ Seurat Processing

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
setwd("~/Dropbox/10x_OSKM_2/3.Seurat")
dir.create("plots_09032020")
dir.create("data_09032020")

#======================================================================================
# Create Seurat object and basic processing
#======================================================================================
# Mapped to custom genome (mm10, FUW-M2rtTA full-length, and each 2A sequence)

#=======================================================================
#=======================================================================
# 10x Lane 1 ----
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data1 <- Read10X("../1.CellRanger/output1_09022020") #31057 x 3006, genes x cells

# Read in LMO count matrix, remove unmapped row
lmo.data1 <- Read10X("../2.LMO/output1_09032020/umi_count/", gene.column=1) #5 x 3006, samples and unmapped x cells
lmo.data1 <- lmo.data1[1:4,] #4 x 3006, samples x cells

# Remove any cells without barcodes and barcodes without cells (none)
colnames(svz.data1) <- sub("-1","",colnames(svz.data1))
barcode.intersect <- intersect(colnames(svz.data1), colnames(lmo.data1)) #3006
lmo.data1 <- lmo.data1[ , barcode.intersect] #4 x 3006
svz.data1 <- svz.data1[ , barcode.intersect] #31057 x 3006

# Make Seurat Object
svz1 <- CreateSeuratObject(counts = svz.data1, project = "Lane1") #31057 features across 3006 samples within 1 assay 
svz1 <- NormalizeData(svz1)
svz1 <- FindVariableFeatures(svz1)
svz1 <- ScaleData(svz1)


#======================================================================================
# QC Metrics and Filtering
#======================================================================================

svz1[["percent.mt"]] <- PercentageFeatureSet(svz1, pattern = "^mt-")

# cut off of 1000 UMI instead of elbow from cell ranger; less than 10% mito; more than 500 features.
svz1 <- subset(svz1, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 
#31061 features across 2760 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots_09032020/violin.mt.feature.counts.fil.lane1.pdf", width = 16, height = 10)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 17 PCs for tSNE visualization and graph-based clustering
svz1 <- RunPCA(svz1, verbose = FALSE)
ElbowPlot(svz1, ndims = 50)
ggsave("plots_09032020/elbow.lane1.pdf", width=7.3, height = 4.3)

svz1 <- FindNeighbors(svz1, dims = 1:23)
svz1 <- FindClusters(svz1, resolution = 0.15)
svz1 <- RunUMAP(svz1, dims = 1:23)
DimPlot(svz1)
ggsave("plots_04012020/umap.initial.fil.pdf", width = 10, height = 10)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

svz1[["LMO"]] <- CreateAssayObject(counts = lmo.data1)
svz1 <- NormalizeData(object = svz1, assay = "LMO", normalization.method = "CLR")
svz1 <- HTODemux(svz1, assay = "LMO", positive.quantile = 0.93)
table(svz1$LMO_classification.global)
# Doublet Negative  Singlet 
# 547      821     1638  

color_pal.6 <- tableau_color_pal(palette = "Tableau 10")(6)

# Group cells based on the max HTO signal
Idents(svz1) <- "LMO_maxID"
RidgePlot(svz1, assay = "LMO", features = rownames(svz1[["LMO"]])[1:4], ncol = 2)
ggsave("plots_09032020/ridge.lmo.fil.CLR.lane1.pdf", height = 11, width = 25)

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz1, assay = "LMO", ncells = 2000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots_09032020/heatmap.lmo.lane1.png", width = 10.9, height=2.88)


#=======================================================================
#=======================================================================
# 10x Lane 2
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data2 <- Read10X("../1.CellRanger/output2_09022020") #31057 x 2715, genes x cells

# Read in LMO count matrix, remove unmapped row
lmo.data2 <- Read10X("../2.LMO/output2_09032020/umi_count/", gene.column=1) #5 x 2715, samples and unmapped x cells
lmo.data2 <- lmo.data2[1:4,] #4 x 2715, samples x cells

# Remove any cells without barcodes and barcodes without cells
colnames(svz.data2) <- sub("-1","",colnames(svz.data2))
barcode.intersect <- intersect(colnames(svz.data2), colnames(lmo.data2)) #2715
lmo.data2 <- lmo.data2[ , barcode.intersect] #4 x 2715
svz.data2 <- svz.data2[ , barcode.intersect] #31057 x 2715

# Make Seurat object
svz2 <- CreateSeuratObject(counts = svz.data2, project = "Lane2") #31057 features across 2715 samples within 1 assay 
svz2 <- NormalizeData(svz2)
svz2 <- FindVariableFeatures(svz2)
svz2 <- ScaleData(svz2)


#======================================================================================
# QC Metrics and Filtering
#======================================================================================

svz2[["percent.mt"]] <- PercentageFeatureSet(svz2, pattern = "^mt-")

# cut off of 1000 UMI instead of elbow from cell ranger; less than 10% mito; more than 500 features.
svz2 <- subset(svz2, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 
#31057 features across 2528 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots_09032020/violin.mt.feature.counts.fil.lane2.pdf", width = 16, height = 10)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 23 PCs for tSNE visualization and graph-based clustering
svz2 <- RunPCA(svz2, verbose = FALSE)
ElbowPlot(svz2, ndims = 50)
ggsave("plots_09032020/elbow.lane2.pdf", width=7.3, height = 4.3)

svz2 <- FindNeighbors(svz2, dims = 1:23)
svz2 <- FindClusters(svz2, resolution = 0.15)
svz2 <- RunUMAP(svz2, dims = 1:23)
DimPlot(svz2)
ggsave("plots_09032020/umap.initial.fil.lane2.pdf", width = 10, height = 10)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz2), colnames(lmo.data2)) #2528
lmo.data2 <- lmo.data2[ , barcode.intersect] #4 x 2528
svz2[["LMO"]] <- CreateAssayObject(counts = lmo.data2)
svz2 <- HTODemux(svz2, assay = "LMO", positive.quantile = 0.93)
# Doublet Negative  Singlet 
# 305      700     1523 

color_pal.6 <- tableau_color_pal(palette = "Tableau 10")(6)

# Group cells based on the max HTO signal
Idents(svz2) <- "LMO_maxID"
RidgePlot(svz2, assay = "LMO", features = rownames(svz2[["LMO"]])[1:4], ncol = 2)
ggsave("plots_09032020/ridge.lmo.fil.CLR.lane2.pdf", height = 11, width = 25)

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz2, assay = "LMO", ncells = 5000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots_09032020/heatmap.lmo.lane2.png", width = 10.9, height=2.88)



#=======================================================================
#=======================================================================
# 10x Lane 3
#=======================================================================
#=======================================================================
# Read in Cellranger count matrix and create basic seurat object
svz.data3 <- Read10X("../1.CellRanger/output3_09022020") #31057 x 8182, genes x cells

# Read in LMO count matrix, remove unmapped row
lmo.data3 <- Read10X("../2.LMO/output3_09032020/umi_count/", gene.column=1) #5 x 8182, samples and unmapped x cells
lmo.data3 <- lmo.data3[1:4,] #4 x 8182, samples x cells

# Remove cells without barcodes and barcodes without cells
colnames(svz.data3) <- sub("-1","",colnames(svz.data3))
barcode.intersect <- intersect(colnames(svz.data3), colnames(lmo.data3)) #8182
lmo.data3 <- lmo.data3[ , barcode.intersect] #4 x 8182
svz.data3 <- svz.data3[ , barcode.intersect] #31057 x 8182

# Make Seurat object
svz3 <- CreateSeuratObject(counts = svz.data3, project = "Lane3") #31057 features across 8182 samples within 1 assay 
svz3 <- NormalizeData(svz3)
svz3 <- FindVariableFeatures(svz3)
svz3 <- ScaleData(svz3)


#======================================================================================
# QC Metrics and Filtering
#======================================================================================

svz3[["percent.mt"]] <- PercentageFeatureSet(svz3, pattern = "^mt-")

# cut off of 1000 UMI instead of elbow from cell ranger; less than 10% mito; more than 500 features.
svz3 <- subset(svz3, subset = nCount_RNA > 1000 & percent.mt < 10 & nFeature_RNA > 500) 
#31057 features across 7743 samples 

# Visualize QC metrics as a violin plot
VlnPlot(svz3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.4)
ggsave("plots_09032020/violin.mt.feature.counts.fil.lane3.pdf", width = 16, height = 10)


#======================================================================================
# Dim reduction and clustering
#======================================================================================

# Run PCA, select 23 PCs for tSNE visualization and graph-based clustering
svz3 <- RunPCA(svz3, verbose = FALSE)
ElbowPlot(svz3, ndims = 50)
ggsave("plots_09032020/elbow.lane3.pdf", width=7.3, height = 4.3)

svz3 <- FindNeighbors(svz3, dims = 1:23)
svz3 <- FindClusters(svz3, resolution = 0.15)
svz3 <- RunUMAP(svz3, dims = 1:23)
DimPlot(svz3)
ggsave("plots_09032020/umap.initial.fil.lane3.pdf", width = 10, height = 10)


#======================================================================================
# Add LMO sample label data to object
#======================================================================================

barcode.intersect <- intersect(colnames(svz3), colnames(lmo.data3)) #7743
lmo.data3 <- lmo.data3[ , barcode.intersect] #4 x 7743
svz3[["LMO"]] <- CreateAssayObject(counts = lmo.data3)

svz3 <- NormalizeData(object = svz3, assay = "LMO", normalization.method = "CLR")
svz3 <- HTODemux(svz3, assay = "LMO", positive.quantile = 0.94)
# Doublet Negative  Singlet 
# 1659     3026     3058  

color_pal.6 <- tableau_color_pal(palette = "Tableau 10")(6)

# Group cells based on the max HTO signal
Idents(svz3) <- "LMO_maxID"
RidgePlot(svz3, assay = "LMO", features = rownames(svz3[["LMO"]])[1:4], ncol = 2)
ggsave("plots_09032020/ridge.lmo.fil.CLR.lane3.pdf", height = 11, width = 25)

# Heatmap version of above ridgeplots (simple but less quant)
HTOHeatmap(svz3, assay = "LMO", ncells = 5000) + scale_fill_viridis(option="viridis", direction=1)
ggsave("plots_09032020/heatmap.lmo.lane3.png", width = 10.9, height=2.88)



#======================================================================================
# Save Seurat objects
#======================================================================================
saveRDS(svz1, paste0("data_09032020/svz1.", Sys.Date(), ".rds"))
saveRDS(svz2, paste0("data_09032020/svz2.", Sys.Date(), ".rds"))
saveRDS(svz3, paste0("data_09032020/svz3.", Sys.Date(), ".rds"))

svz1 <- readRDS("data_09032020/svz1.2020-09-07.rds")
svz2 <- readRDS("data_09032020/svz2.2020-09-07.rds")
svz3 <- readRDS("data_09032020/svz3.2020-09-07.rds")


#======================================================================================
# Merge Seurat objects
#======================================================================================

svz <- merge(svz1, c(svz2, svz3), add.cell.ids = c("l1","l2","l3"), project = "OSKM2", merge.data = TRUE)
#31069 features across 13031 samples within 2 assays 

svz <- NormalizeData(svz)
svz <- FindVariableFeatures(svz)
svz <- ScaleData(svz)

# Run PCA
svz <- RunPCA(svz, verbose = FALSE)
ElbowPlot(svz, ndims = 50)
ggsave("plots_09032020/elbow.combined.pdf", width=7.3, height = 4.3)

svz <- FindNeighbors(svz, dims = 1:21)
svz <- FindClusters(svz, resolution = 0.15)
svz <- RunUMAP(svz, dims = 1:21)
DimPlot(svz)
ggsave("plots_09032020/umap.initial.fil.combined.pdf", width = 10, height = 10)




#======================================================================================
# Subset to singlets
#======================================================================================

# Set colors
color_pal.1 <- tableau_color_pal(palette = "Tableau 10")(1)
color_pal.4 <- tableau_color_pal(palette = "Tableau 10")(4)

# Remove doublets and Negatives
svz_sing <- subset(svz, subset = LMO_classification.global == "Singlet") 
#31069 features across 6097 samples within 2 assays 

# UMAP based on transcriptome reduction, colored by doublet/singlet status
DimPlot(svz_sing, reduction = "umap")
ggsave("plots_09032020/umap.sing.combined.pdf", width = 5.9, height = 4.6)

# UMAP based on transcriptome reduction, colored by sample barcode
DimPlot(svz_sing, reduction = "umap", group.by="hash.ID")
ggsave("plots_09032020/umap.sing.combined.sample.pdf", width = 9.2, height = 4.6)

# UMAP based on transcriptome reduction, split by lane
DimPlot(svz_sing, reduction = "umap", split.by="orig.ident")
ggsave("plots_09032020/umap.sing.combined.lanesplit.pdf", width = 9.2, height = 4.6)




#======================================================================================
# Rerun processing on subset
#======================================================================================

# Normalize, standardize.  
svz_sing <- SCTransform(svz_sing)

# Reduce
svz_sing <- RunPCA(svz_sing, verbose = FALSE)
ElbowPlot(svz_sing, ndims = 50)
ggsave("plots_09032020/elbow.sing.combined.sct.pdf", width = 5.3, height = 2.8)

svz_sing <- FindNeighbors(svz_sing, dims = 1:20)
svz_sing <- FindClusters(svz_sing, resolution = 0.5)
svz_sing <- RunUMAP(svz_sing, dims = 1:20)

DimPlot(svz_sing, pt.size = .6)
ggsave("plots_09032020/umap.sing.combined.sct.pdf", width = 5.9, height = 4.6)

DimPlot(svz_sing, group.by = "hash.ID", pt.size = .6)
ggsave("plots_09032020/umap.sing.combined.sct.sample.pdf", width = 10, height = 6)

DimPlot(svz_sing, group.by = "orig.ident", pt.size = .6)
ggsave("plots_09032020/umap.sing.combined.sct.lane.pdf", width = 10, height = 6)

DimPlot(svz_sing, split.by = "orig.ident", pt.size = .6)
ggsave("plots_09032020/umap.sing.combined.sct.lanesplit.pdf", width = 10, height = 6)


#======================================================================================
# Save Seurat object
#======================================================================================
saveRDS(svz_sing, paste0("data_09032020/svz_sing.combined.", Sys.Date(), ".rds"))
svz_sing <- readRDS("data_09032020/svz_sing.combined.2020-09-08.rds")


#======================================================================================
# Find markers
#======================================================================================

# Find marker genes for clusters and save
svz_sing.markers <- FindAllMarkers(object=svz_sing)
write.csv(svz_sing.markers, paste0('data_09032020/svz_sing.combined.markers.', Sys.Date(), '.csv'))
saveRDS(svz_sing.markers, paste0("data_09032020/svz_sing.combined.markers_", Sys.Date(), ".rds"))

# Pull out top 10 for heatmap 
top10 <- svz_sing.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz_sing, features = top10$gene, lines.width = 10, label = F) + NoLegend() + theme(axis.text.y = element_text(size = 10))
ggsave("plots_09032020/heatmap.markers.top10.png", height = 20, width = 20)
write.csv(top10, paste0('data_09032020/top10_svz_sing.combined.markers.', Sys.Date(), '.csv'))


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