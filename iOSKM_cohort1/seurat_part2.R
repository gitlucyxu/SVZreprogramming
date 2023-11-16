# Whole-body (iOSKM) partial reprogramming, cohort 1
# LMO + 6 SVZ Seurat Processing

# PART TWO

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
library(RColorBrewer)
sessionInfo()


setwd("~/Dropbox/OSKM_10x/3.Seurat") # LOCAL
#dir.create("plots2_04012020")

#======================================================================================
# Load data
markers <- readRDS("data_04012020/svz_sing.0.8.markers_2020-04-12.rds")
svz_sing_0.8 <- readRDS("data_04012020/svz_sing.0.8_2020-04-12.rds")
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)
t4 <- tableau_color_pal(palette = "Tableau 10")(4)[c(2,3,4,1)]

svz_sing_0.8 <- RunUMAP(svz_sing_0.8, dims = 1:18, min.dist = .8, spread = .5, seed.use = 221)
DimPlot(svz_sing_0.8, group.by = "seurat_clusters", cols = t20, label = TRUE, label.size = 5, pt.size = .6)
ggsave(paste0("plots2_07012020/umap.seed", 221, ".pdf"), width=9, height=6)
colnames(svz_sing_0.8[[]])


#======================================================================================
# Assign celltypes
#======================================================================================

# Increase resolution to separate small but distinct clusters
svz_sing <- FindClusters(svz_sing_0.8, resolution = 5.6)
svz_sing <- RunUMAP(svz_sing, dims = 1:18, min.dist = .8, spread = .5, seed.use = 221)
DimPlot(svz_sing, group.by = "seurat_clusters", label = TRUE, label.size = 5, pt.size = .6)
ggsave(paste0("plots2_07012020/umap.38clusters.seed", 221, ".pdf"), width=9, height=6)

# Recombine large clusters
clusternum <- plyr::mapvalues(x = svz_sing$seurat_clusters, 
                             from = c("0","1","2","4","16","20","31","6","26","19","10","28",
                                      "22","13","11","7","29","3","24","18","23", "8",
                                      "25","21","9",
                                      "14","5",
                                      "15","27"
                                      ), 
                             to =   c("A","A","A","A","A","A","A","A","A","A","A","A",
                                      "B","B","B","B","B","B","B","B","B","B",
                                      "C","C","C",
                                      "D","D",
                                      "E","E"
                                      )
)

svz_sing$NewClusters <- clusternum
#Levels: A B D C 12 E 17 30 32 33 34 35 36 37

new.cluster.ids <- c("Microglia",
                     "Oligodendrocyte",
                     "Astrocyte_qNSC",
                     "Endothelial",
                     "qNSC.primed",
                     "Neuroblast",
                     "aNSC_NPC",
                     "30",
                     "32",
                     "33",
                     "34",
                     "35",
                     "36",
                     "37"
                     )

Idents(svz_sing) <- svz_sing$NewClusters
names(new.cluster.ids) <- levels(svz_sing)
svz_sing <- RenameIdents(svz_sing, new.cluster.ids)
svz_sing[["Celltype.38collapsed"]] <- Idents(svz_sing)
unique(svz_sing@meta.data$Celltype.38collapsed) 
#Levels: Microglia Oligodendrocyte Astrocyte_qNSC Endothelial qNSC.primed Neuroblast aNSC_NPC 30 32 33 34 35 36 37

DimPlot(svz_sing, group.by = "Celltype.38collapsed", cols = t20, label = TRUE, label.size = 5, pt.size = .6)
ggsave("plots2_07012020/umap.celltype.38collapsed.pdf", width=9, height=6)

# Find markers
newmarkers <- FindAllMarkers(object=svz_sing)
write.csv(newmarkers, paste0('data_07012020/svz_sing.38clusterscollapsed.markers_', Sys.Date(), '.csv'))
saveRDS(newmarkers, paste0("data_07012020/svz_sing.38clusterscollapsed.markers_", Sys.Date(), ".rds"))
newmarkers <- readRDS("data_07012020/svz_sing.38clusterscollapsed.markers_2020-06-29.rds")

top10 <- newmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(svz_sing, features = top10$gene, group.by = "Celltype.38collapsed", group.colors = t20, lines.width = 3, label = F,
          cells = WhichCells(svz_sing, downsample = 50)) +
  NoLegend() + theme(axis.text.y = element_text(size = 10))

DoHeatmap(svz_sing, features = top10$gene, group.by = "Celltype.38collapsed", group.colors = t20, lines.width = 2, label = F,
          cells = WhichCells(svz_sing, idents = c("30","32","33","34","35","36","37"))) +
          NoLegend() + theme(axis.text.y = element_text(size = 10))

write.csv(top10, paste0('data_07012020/top10_svz_sing.38clusterscollapsed.markers_', Sys.Date(), '.csv'))
saveRDS(svz_sing, paste0("data_07012020/svz_38clusterscollapsed_celltypes_", Sys.Date(), ".rds"))

# 34 + 35: Mural. 34 are vascular smooth muscle cells (Acta2, Myl9, Pdlim3) and 35 are pericytes (Kcnj8, Ifitm1,Rgs5)
# Other marker genes from https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
FeaturePlot(svz_sing, features = c("Pdlim3","Acta2","Myl9","Kcnj8","Rgs5" ,"Ifitm1","Pdgfrb"))
ggsave("plots2_07012020/featureplot.mural.pdf", height = 9, width = 10)

# 30: OPCs. some Pdgfra+. Opalin-Mog-.
FeaturePlot(svz_sing, features = c("Pdgfra","Opalin","Mbp","Mog"))
ggsave("plots2_07012020/featureplot.OPC.pdf", height = 6, width = 6)

# 32, 33: Microglia/macrophage subtypes
# Cx3cr1 is a marker for interstitial microglia, Mrc1 is a marker for perivascular microglia and macrophages. https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
# 32 is Cx3cr1-Mrc1+. Assigning as macrophage
# 33 is Cx3cr1-Mrc1-. Compared to cluster 32, Ifi27l2a low, Lyz2 low.
FeaturePlot(svz_sing, features = c("Cx3cr1", "Mrc1","Ifi27l2a", "Lyz2"))
ggsave("plots2_07012020/featureplot.microglia.pdf", height = 6, width = 6)

FeaturePlot(svz_sing, features = c("Rpl26", "Rplp1", "Malat1"), ncol = 3)
ggsave("plots2_07012020/featureplot.cluster33.pdf", height = 3, width = 10)

# 32: Monocyte/macrophage-specific markers from https://doi-org.stanford.idm.oclc.org/10.1016/j.immuni.2018.11.004 
FeaturePlot(svz_sing,features=c("F13a1","H2-Aa","Ccr2","Lyve1","Mgl2", "P2ry12"))
ggsave("plots2_07012020/featureplot.macrophage-Stevens.pdf", height = 9, width = 6)

# 36: unknown. enriched for proliferation genes. 
FeaturePlot(svz_sing, features = c("Birc5", "Cenpa","P2ry12","Ifi27l2a"))
ggsave("plots2_07012020/featureplot.cluster36.pdf", height = 6, width = 6)

# 37: Ependymal cells
# aSMA (acta2) is used as a marker to sort. Other ependymal markers are also specific to this cluster. https://doi-org.stanford.idm.oclc.org/10.1016/j.cell.2018.03.063 
FeaturePlot(svz_sing, features = c("Acta2", "Rarres2", "Tmem212", "Mia"))
ggsave("plots2_07012020/featureplot.ependymal.pdf", height = 6, width = 6)


#---------------------------
# Add new cell type annotations 
unique(svz_sing@meta.data$Celltype.38collapsed)
celltypes <- plyr::mapvalues(x = svz_sing$Celltype.38collapsed, 
                             from = c("30","32","33","34","35","36","37"), 
                             to = c("OPC", "Macrophage", "33", "Vasc.Smooth.Muscle", "Pericyte", "36","Ependymal"))
svz_sing$Celltype.38collapsed <- celltypes
unique(svz_sing@meta.data$Celltype.38collapsed)

# Low res option: combines mural cells together, qNSC-primed in with qNSCs, and unassigned clusters 33 and 36 with microglia.
svz_sing$Celltype.LowRes <- svz_sing$Celltype.38collapsed
lowres <- plyr::mapvalues(x = svz_sing$Celltype.38collapsed, 
                             from = c("qNSC.primed","33","36","Vasc.Smooth.Muscle","Pericyte"), 
                             to = c("Astrocyte_qNSC", "Microglia", "Microglia", "Mural", "Mural"))
svz_sing$Celltype.LowRes <- lowres

Idents(svz_sing) <- svz_sing$Celltype.38collapsed
svz_sing[["Celltype"]] <- Idents(svz_sing)
Idents(svz_sing) <- svz_sing@meta.data$Celltype


# DEFINE COLORS BY CELLTYPE --------------------
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal")
svz_sing$Celltype.LowRes <- factor(svz_sing$Celltype.LowRes,  levels=CELLS, ordered=T)

new_colors <- c("#03c03c", "#4E79A7", "#966fd6",
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#db7093", "#2F4F4F")
names(new_colors) <- levels(svz_sing$Celltype.LowRes)
# ----------------------------------------

DimPlot(svz_sing, group.by = "Celltype.LowRes", label = TRUE, cols = alpha(new_colors,0.7), label.size = 4, pt.size = .6, repel = T)
ggsave(paste0("plots2_07012020/umap.celltypes.lowres.seed", 221, ".labels.pdf"), width=9, height=6)

DimPlot(svz_sing, group.by = "Celltype.LowRes", cols = alpha(new_colors,0.7), pt.size = .6)
ggsave(paste0("plots2_07012020/umap.celltypes.lowres.seed", 221, ".pdf"), width=9, height=6)


#====================================
# Add treatment info
#====================================
treatment <- plyr::mapvalues(x = svz_sing$hash.ID, 
                        from = c("Sample1-961L-DOB20180508-2Dox0-CTCTAGAC", "Sample2-955L-DOB20180501-untr-ACCAATGC", "Sample3-928L-DOB20180327-2Dox5-AGTTGCGT", "Sample4-960L-DOB20180508-2Dox0-CGAACAAG", "Sample5-929L-DOB20180327-2Dox5-GTACCTGT", "Sample6-953L-DOB20180501-untr-GAAGCTTG"), 
                        to = c("2Dox0", "untr", "2Dox5", "2Dox0", "2Dox5", "untr")
)
svz_sing$Treatment <- treatment

# Adjust factors
TREATMENT <- c("2Dox0", "2Dox5", "untr")
svz_sing$Treatment <- factor(svz_sing$Treatment,  levels=TREATMENT, ordered=T)


#======================================================================================
# SAVE
#======================================================================================
saveRDS(svz_sing, paste0("data_07012020/svz_celltypes_", Sys.Date(), ".rds"))


#=================
# UMAPs 
#=================
setwd("~/Dropbox/OSKM_10x/3.Seurat") # LOCAL

col3 <- c("turquoise2","firebrick", "mediumorchid", "darkslategray3")
names(col3) <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")

svz_sing <- readRDS(paste0("data_07012020/svz_celltypes_2020-06-29.rds"))

svz_sing$Age_Treatment <- plyr::mapvalues(svz_sing$Treatment, 
                                          from = c("untr", "2Dox0", "2Dox5"),
                                          to = c("old_untr","old_2Dox0","old_2Dox5"))

svz_sing$Age_Treatment <- factor(svz_sing$Age_Treatment, levels = c("old_untr","old_2Dox0","old_2Dox5"), ordered = T)

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal")
svz_sing$Celltype <- factor(svz_sing$Celltype.LowRes, levels = CELLS, ordered = T)

newcolors <- c("#03c03c", "#4E79A7", "#966fd6",
               "#ffdf00", "#B6992D", "#FFBE7D", 
               "#aec6cf", "#86BCB6", "#e5aa70", "#2F4F4F")
newcolors <- alpha(newcolors, 0.7)
names(newcolors) <- levels(svz_sing$Celltype)

table(svz_sing$Age_Treatment)
# old_untr old_2Dox0 old_2Dox5 
# 550       644      1056


# NB: The samples labled 'old_2Dox5' are not included in this paper. 

png('plots2_07012020/forpaper_umap.no2dox5.celltypes.nolegend.png', width=3, height=3, units = 'in', res=300)
DimPlot(svz_sing %>% subset(Age_Treatment!="old_2Dox5"), reduction = "umap", group.by = "Celltype", label = F, pt.size = .8, cols = newcolors, shape.by = "Treatment") + NoLegend() +
  scale_shape_manual(values=c(16,16,16), guide="none") +
  theme(axis.title = element_text(size=12)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 12, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 12, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  scale_x_continuous(breaks=seq(-10,10,10)) +
  theme(legend.text = element_text(size=10.5))
dev.off()

png('plots2_07012020/forpaper_umap.no2dox5.treatment.nolegend.png', width=3, height=3, units = 'in', res=400)
DimPlot(svz_sing %>% subset(Age_Treatment!="old_2Dox5"), reduction = "umap", group.by = "Age_Treatment", label = F, pt.size = .8, cols = alpha(col3,0.6), shape.by = "Age_Treatment") + NoLegend() +
  scale_shape_manual(values=c(16,16,16), guide="none") +
  theme(axis.title = element_text(size=12)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 12, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 12, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  scale_x_continuous(breaks=seq(-10,10,10)) 
dev.off()


# Downsample ----
Idents(svz_sing) <- "Treatment"
downsampled <- subset(svz_sing, downsample = 550)
table(downsampled$Treatment)

downsampled$Age_Treatment <- factor(downsampled$Age_Treatment, levels = c("old_untr","old_2Dox0","old_2Dox5"), ordered = T)

p <- DimPlot(downsampled %>% subset(Age_Treatment!="old_2Dox5"), reduction = "umap", group.by = "Age_Treatment", split.by = "Age_Treatment", label = F, pt.size = .8, cols = alpha(col3,0.6), shape.by = "Age_Treatment") + NoLegend() +
  scale_shape_manual(values=c(16,16,16), guide="none") +
  theme(axis.title = element_text(size=12)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 12, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 12, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  theme(strip.text = element_blank()) +
  scale_x_continuous(breaks=seq(-10,10,10)) +
  theme(legend.text = element_text(size=10.5))

png('plots2_07012020/forpaper_umap.no2dox5.treatment.split.png', width=4, height=3, units = 'in', res=300)
p
dev.off()




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
