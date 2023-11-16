# Whole-body (iOSKM) partial reprogramming, cohort 2
# LMO + 12 SVZ Seurat Processing
# PART TWO

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

setwd("~/Dropbox/10x_OSKM_2/3.Seurat") # LOCAL
dir.create("plots2_09032020")

#======================================================================================
# Load data
markers <- readRDS("data_09032020/svz_sing.combined.markers_2020-09-08.rds")
svz_sing <- readRDS("data_09032020/svz_sing.combined.2020-09-08.rds")
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)
t4 <- tableau_color_pal(palette = "Tableau 10")(4)[c(2,3,4,1)]

svz_sing <- RunUMAP(svz_sing, dims = 1:20)
DimPlot(svz_sing, group.by = "seurat_clusters", label = TRUE, label.size = 5, pt.size = .6)
ggsave("plots_09032020/umap.sing.combined.sct.labels.pdf", width=7, height=6)


#====================================
# Assign celltypes
#====================================

# Increase cluster resolution to separate small clusters
svz_sing <- FindNeighbors(svz_sing, dims = 1:20)
svz_sing <- FindClusters(svz_sing, resolution = 1.5)
svz_sing <- RunUMAP(svz_sing, dims = 1:20, min.dist = .5, spread = .75, seed.use = 2)
DimPlot(svz_sing, group.by = "seurat_clusters", label = TRUE, label.size = 5, pt.size = .6)
ggsave(paste0("plots2_09032020/umap.30clusters.seed", 22, ".pdf"), width=7, height=6)

# Find marker genes for clusters and save
svz_sing.markers <- FindAllMarkers(object=svz_sing)
write.csv(svz_sing.markers, paste0('data_09032020/svz_sing.combined.31clusters.markers.', Sys.Date(), '.csv'))
saveRDS(svz_sing.markers, paste0("data_09032020/svz_sing.combined.31clusters.markers_", Sys.Date(), ".rds"))

# Pull out top 10 for heatmap 
top10 <- svz_sing.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(svz_sing, features = top10$gene, lines.width = 10, label = T) + NoLegend() + theme(axis.text.y = element_text(size = 10))
ggsave("plots2_09032020/heatmap.markers.31clusters.top10.png", height = 20, width = 20)
write.csv(top10, paste0('data_09032020/top10_svz_sing.combined.31clusters.markers.', Sys.Date(), '.csv'))

#--------
# Marker expression
#--------

# NSC lineage 
genes <- c("Clu", "Aldoc", "Meg3", "Tubb2b", "Pclaf", "Cenpf")
FeaturePlot(svz_sing, features = genes)
ggsave("plots2_09032020/umap.NSClineage-markers.pdf", width = 8, height = 10)

# qNSC clusters 0,13,18 
genes <- c("Aldoc", "Clu","Plpp3","Mt3","Slc1a2")
genes <-c("Thbs4","Fxyd1","Cldn10","Aldoc","Ntsr2","Cspg5","Mlc1","Clu","Bcan","Chchd10")
FeaturePlot(svz_sing,features=genes, order = T)

# cluster 15 aNSCs
genes <-c("Egfr","Ascl1","Mcm2","Cenpm","Cdc6","E2f1")
ggsave("plots2_09032020/featureplot.markers.15.png", height = 10, width = 10)

#cluster 8 aNSCs
genes <-c("Hist1h1b","Hist1h2ae","Hist1h3c","Top2a","Hist1h2ap","Pclaf","Hist1h2ab","Mki67","H2afx","Cenpf")

# Microglia, olig, and endothelial marker genes
genes <- c("C1qc", "Junb", "Cd52", "Opalin", "Tspan2", "Klk6", "Cxcl12", "Ly6c1", "Itm2a")
FeaturePlot(svz_sing, features = genes)
ggsave("plots2_09032020/umap.micro-olig-endo-markers.pdf", width = 12, height = 10)

# Mural cells
# cluster 26: vascular smooth muscle cells (Acta2, Tagln, Myh11). 
# cluster 20: pericytes (Kcnj8,Ifitm1,Vtn)
# Other marker genes from https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
FeaturePlot(svz_sing, features = c("Acta2","Tagln","Myh11","Vtn","Kcnj8","Ifitm1","Pdgfrb"))
ggsave("plots2_09032020/featureplot.markers.26vsmc.20pericyte.mural.pdf", height = 9, width = 10)

# cluster 22: OPCs. Pdgfra+, Opalin-Mog-.
FeaturePlot(svz_sing, features = c("Pdgfra","Opalin","Mbp","Mog"), order= T)
ggsave("plots2_09032020/featureplot.markers.22OPC.pdf", height = 6, width = 6)

# Microglia/macrophage subtypes
# cluster 24: macrophage (Cx3cr1-Mrc1+)
# Cx3cr1 is a marker for interstitial microglia, Mrc1 is a marker for perivascular microglia and macrophages. https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
FeaturePlot(svz_sing, features = c("Cx3cr1", "Mrc1","Ifi27l2a", "Lyz2"))
ggsave("plots2_09032020/featureplot.markers.24macrophage.pdf", height = 6, width = 6)

# Monocyte/macrophage-specific markers from https://doi-org.stanford.idm.oclc.org/10.1016/j.immuni.2018.11.004 
FeaturePlot(svz_sing,features=c("F13a1","H2-Aa","Ccr2","Lyve1","Mgl2", "P2ry12"), order = T)
ggsave("plots2_09032020/featureplot.markers.macrophage-Stevens.pdf", height = 9, width = 6)

# cluster 25: Ependymal cells
# aSMA (acta2) is used as a marker to sort. Other ependymal markers are also specific to this cluster. https://doi-org.stanford.idm.oclc.org/10.1016/j.cell.2018.03.063. Top GO term is cilium movement. 
FeaturePlot(svz_sing, features = c("Acta2", "Rarres2", "Tmem212", "Mia","Dynlrb2","Ccdc153"), order = T)
ggsave("plots2_09032020/featureplot.markers.25ependymal.pdf", height = 9, width = 6)

# T cell genes: Cd3, Cd8, Cd4, interferon, activation, and tissue retention genes
genes <- c("Cd3g", "Cd3d", "Cd3e", "Cd8a", "Cd4", "Ifng", "Cd69", "Xcl1", "Itga4") #"Cd8b"
FeaturePlot(svz_sing, features = genes, order = TRUE)
ggsave("plots2_09032020/umap.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)


#--------
# Assign celltypes
#--------
#31 clusters with resolution = 1.5
new.cluster.ids <- c("Astrocyte_qNSC_1",
                     "Oligodendrocyte_1",
                     "Microglia_1",
                     "Microglia_2",
                     "Microglia_3",
                     "Oligodendrocyte_2",
                     "Oligodendrocyte_3",
                     "Neuroblast_1",
                     "aNSC_NPC_1",
                     "Endothelial_1",
                     "Microglia_4",
                     "Oligodendrocyte_4",
                     "Neuroblast_2",
                     "Astrocyte_qNSC_2",
                     "Microglia_5",
                     "aNSC_NPC_2",
                     "Neuroblast_3",
                     "Oligodendrocyte_5",
                     "Astrocyte_qNSC_3",
                     "Endothelial_2",
                     "Pericyte",
                     "aNSC_NPC_3",
                     "OPC",
                     "Microglia_6",
                     "Macrophage",
                     "Ependymal",
                     "Vasc.Smooth.Musc",
                     "27",
                     "28",
                     "29",
                     "30"
                     )

Idents(svz_sing) <- "seurat_clusters"
names(new.cluster.ids) <- levels(svz_sing)
svz_sing <- RenameIdents(svz_sing, new.cluster.ids)
svz_sing[["Celltype.31"]] <- Idents(svz_sing)
levels(unique(svz_sing@meta.data$Celltype.31))

# Assign lower resolution celltypes
low.res.ids <-  c("Astrocyte_qNSC",
                    "Oligodendrocyte",
                    "Microglia",
                    "Microglia",
                    "Microglia",
                    "Oligodendrocyte",
                    "Oligodendrocyte",
                    "Neuroblast",
                    "aNSC_NPC",
                    "Endothelial",
                    "Microglia",
                    "Oligodendrocyte",
                    "Neuroblast",
                    "Astrocyte_qNSC",
                    "Microglia",
                    "aNSC_NPC",
                    "Neuroblast",
                    "Oligodendrocyte",
                    "Astrocyte_qNSC",
                    "Endothelial",
                    "Mural",
                    "aNSC_NPC",
                    "OPC",
                    "Microglia",
                    "Macrophage",
                    "Ependymal",
                    "Mural",
                    "27",
                    "28",
                    "29",
                    "30"
)

Idents(svz_sing) <- svz_sing@meta.data$Celltype.31 
names(low.res.ids) <- levels(svz_sing)
svz_sing <- RenameIdents(svz_sing, low.res.ids)
svz_sing[["Celltype.LowRes"]] <- Idents(svz_sing)
unique(svz_sing@meta.data$Celltype.LowRes)
Idents(svz_sing) <- svz_sing@meta.data$Celltype.LowRes
#Levels: Astrocyte_qNSC Oligodendrocyte Microglia Neuroblast aNSC_NPC Endothelial Mural OPC Macrophage Ependymal 27 28 29 30

# Fix factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
              "Oligodendrocyte", "OPC", "Endothelial",
              "Microglia", "Macrophage", "Mural", "Ependymal", "27", "28", "29", "30")
svz_sing$Celltype <- factor(svz_sing$Celltype.LowRes, levels = CELLS, ordered = T)

newcolors <- c("#03c03c", "#4E79A7", "#966fd6",
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#e5aa70", "#2F4F4F", "#59A14F", "#A52A2A","#E15759","#76B7B2")
newcolors <- alpha(newcolors, 0.7)
names(newcolors) <- levels(svz_sing$Celltype)


#======================================================================================
# Celltype markers heatmap 
#======================================================================================

CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Macrophage", "Mural", "Ependymal", "27", "28","29","30")
newcolors <- c("#03c03c", "#4E79A7", "#966fd6",
               "#ffdf00", "#B6992D", "#FFBE7D", 
               "#aec6cf", "#86BCB6", "#e5aa70", "#2F4F4F", "#59A14F", "#A52A2A","#E15759","#76B7B2")
newcolors <- alpha(newcolors, 0.7)
names(newcolors) <- CELLS

# Find marker genes for clusters and save
Idents(svz_sing) <- svz_sing@meta.data$Celltype

svz_sing.markers <- FindAllMarkers(object=svz_sing)
write.csv(svz_sing.markers, paste0('data_09032020/svz_sing.annotated.markers.', Sys.Date(), '.csv'))
saveRDS(svz_sing.markers, paste0("data_09032020/svz_sing.annotated.markers_", Sys.Date(), ".rds"))

# Pull out top 5 for heatmap
top5 <- svz_sing.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Put in order of prevalence
top5$cluster <- factor(top5$cluster, 
                       levels = c("Microglia","Oligodendrocyte","Astrocyte_qNSC","Neuroblast","aNSC_NPC","Endothelial","Mural","OPC","Macrophage","27", "28","29","30"), ordered=T)

Idents(svz_sing) <- factor(svz_sing@meta.data$Celltype, 
                           levels = c("Microglia","Oligodendrocyte","Astrocyte_qNSC","Neuroblast","aNSC_NPC","Endothelial","Mural","OPC","Macrophage","Ependymal","27", "28","29","30"), ordered=T)

DoHeatmap(svz_sing, features = top5$gene, lines.width = 10, label=T, group.colors = newcolors) + NoLegend() + theme(axis.text.y = element_text(size = 20))
ggsave("plots2_09032020/heatmap.markers.annotated.ordered.top5.png", height = 20, width = 20)



#====================================
# Add sample information 
#====================================
unique(svz_sing$hash.ID)
age_genotype <- plyr::mapvalues(x = svz_sing$hash.ID, 
                             from = c("A-L0002-06102018-lane1-2Dox0-TGTGATGG", 
                                      "B-L0003-06102018-lane1-untr-TCAATGGC", 
                                      "C-999L-06102018-lane1-2Dox5-CTCTAGAC", 
                                      "D-L0011-06292018-lane2-untr-ACCAATGC", 
                                      "E-L0004-06292018-lane2-2Dox0-AGTTGCGT", 
                                      "F-BL61-07012018-lane2-untr-CGAACAAG",
                                      "G-BL62-07012018-lane1-2Dox0-GTACCTGT",
                                      "H-BL63-07012018-lane2-2Dox5-GAAGCTTG",
                                      "I-144-04292020-lane3-untr-AAGTACGC",
                                      "J-148-04292020-lane3-untr-ATTCGCAC",
                                      "K-L0034-07072018-lane3-2Dox5-GAGTCGAT",
                                      "L-L0042-08022018-lane3-untr-AAGGCTAG"), 
                             to = c("old_OSKM", 
                                    "old_OSKM", 
                                    "old_OSKM", 
                                    "old_OSKM", 
                                    "old_OSKM", 
                                    "old_BL6",
                                    "old_BL6",
                                    "old_BL6",
                                    "young_OSKM",
                                    "young_OSKM",
                                    "old_OSKM",
                                    "old_OSKM")
)
svz_sing$age_genotype <- age_genotype

age <- plyr::mapvalues(x = svz_sing$age_genotype, 
                                from = c("old_OSKM", 
                                         "old_BL6", 
                                         "young_OSKM"),
                                to = c("old", 
                                       "old", 
                                       "young"))
svz_sing$age <- age

genotype <- plyr::mapvalues(x = svz_sing$age_genotype, 
                       from = c("old_OSKM", 
                                "old_BL6", 
                                "young_OSKM"),
                       to = c("OSKM", 
                              "BL6", 
                              "OSKM"))

svz_sing$genotype <- genotype

treatment <- plyr::mapvalues(x = svz_sing$hash.ID, 
                             from = c("A-L0002-06102018-lane1-2Dox0-TGTGATGG", 
                                      "B-L0003-06102018-lane1-untr-TCAATGGC", 
                                      "C-999L-06102018-lane1-2Dox5-CTCTAGAC", 
                                      "D-L0011-06292018-lane2-untr-ACCAATGC", 
                                      "E-L0004-06292018-lane2-2Dox0-AGTTGCGT", 
                                      "F-BL61-07012018-lane2-untr-CGAACAAG",
                                      "G-BL62-07012018-lane1-2Dox0-GTACCTGT",
                                      "H-BL63-07012018-lane2-2Dox5-GAAGCTTG",
                                      "I-144-04292020-lane3-untr-AAGTACGC",
                                      "J-148-04292020-lane3-untr-ATTCGCAC",
                                      "K-L0034-07072018-lane3-2Dox5-GAGTCGAT",
                                      "L-L0042-08022018-lane3-untr-AAGGCTAG"), 
                             to = c("2Dox0", 
                                    "untr", 
                                    "2Dox5", 
                                    "untr", 
                                    "2Dox0", 
                                    "untr",
                                    "2Dox0",
                                    "2Dox5",
                                    "untr",
                                    "untr",
                                    "2Dox5",
                                    "untr")
)
svz_sing$Treatment <- treatment

# Adjust factors-----------------------------
AGE <- c("young", "old")
svz_sing$age <- factor(svz_sing$age, levels=AGE, ordered=T)
table(svz_sing[[c("hash.ID", "age")]])

GENOTYPE <- c("OSKM", "BL6")
svz_sing$genotype <- factor(svz_sing$genotype, levels=GENOTYPE, ordered=T)
table(svz_sing[[c("hash.ID", "genotype")]])

AGE.GENOTYPE <- c("young_OSKM", "old_OSKM", "old_BL6")
svz_sing$age_genotype <- factor(svz_sing$age_genotype, levels=AGE.GENOTYPE, ordered=T)
table(svz_sing[[c("hash.ID", "age_genotype")]])

TREATMENT <- c("2Dox0", "2Dox5", "untr")
svz_sing$Treatment <- factor(svz_sing$Treatment,  levels=TREATMENT, ordered=T)
table(svz_sing[[c("hash.ID", "Treatment")]])


# SET COLORS -----------------------------
agecolors <- c("gray", "darkred")
names(agecolors) <- levels(svz_sing$age)

agegenotypecolors <- c("gray", "darkred", "black")
names(agegenotypecolors) <- levels(svz_sing$age_genotype)
# ----------------------------------------------------------



#======================================================================================
# SAVE
#======================================================================================

saveRDS(svz_sing, paste0("data_09032020/svz_celltypes_metadata_", Sys.Date(), ".rds"))
svz_sing <- readRDS(paste0("data_09032020/svz_celltypes_metadata_", "2020-09-22", ".rds"))


#======================================================================================
# Subset Seurat objects
#======================================================================================

# Subset samples to mice included in this study, and add metadata
sub.o <- subset(svz_sing, genotype == "OSKM")

age_treatment <- paste0(sub.o@meta.data$age, "_", sub.o@meta.data$Treatment)
m <- data.frame("Age_Treatment" = age_treatment)
rownames(m) <- rownames(sub.o@meta.data)
sub.o <- AddMetaData(object = sub.o, metadata = m)
sub.o$Age_Treatment <- factor(sub.o$Age_Treatment, levels=c("young_untr","old_untr","old_2Dox0","old_2Dox5"), ordered=T)
sub.o$Celltype <- factor(sub.o$Celltype.LowRes, levels = CELLS, ordered = T)

subsub <- subset(sub.o, Age_Treatment!="old_2Dox5") #3857 cells

# Colors
col2 <- c("turquoise2","firebrick", "mediumorchid")
names(col2) <- c("young_untr","old_untr","old_2Dox0")



#======================================================================================
# UMAPs
#======================================================================================

png('plots2_09032020/forpaper_umap.celltypes.sub-OSKM.no2Dox5.nolegend.png', width=4, height=4, units = 'in', res=300)
DimPlot(subsub, reduction = "umap", group.by = "Celltype.LowRes", label = F, pt.size = .6, cols = newcolors, shape.by = "genotype") + NoLegend() +
  scale_shape_manual(values=16, guide="none") +
  theme(axis.title = element_text(size=12)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 12, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 12, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  scale_x_continuous(breaks=seq(-10,10,10)) +
  theme(legend.text = element_text(size=10))
dev.off()


png('plots2_09032020/forpaper_umap.age-treatment.sub-OSKM.no2Dox5.nolegend.png', width=4, height=4, units = 'in', res=300)
DimPlot(subsub, reduction = "umap", group.by = "Age_Treatment", label = F, pt.size = .6, cols = alpha(col2,0.6), shape.by = "genotype") + NoLegend() +
  scale_shape_manual(values=16, guide="none") +
  theme(axis.title = element_text(size=12)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 12, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 12, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  scale_x_continuous(breaks=seq(-10,10,10)) +
  theme(legend.text = element_text(size=10.5))
dev.off()


# Downsample to lowest number of cells per treatment condition

Idents(subsub) <- 'Age_Treatment'
downsampled <- subset(subsub, downsample = 732)

png('plots2_09032020/forpaper_umap.sub-OSKM.no2Dox5.split-age-treatment.downsampled.nolegend.png', width=5.5, height=3, units = 'in', res=300)
DimPlot(downsampled, reduction = "umap", label = F, label.size = 5, pt.size = 0.8, split.by = "Age_Treatment", group.by = "Age_Treatment", shape.by = "genotype", cols = alpha(col3,0.4)) + NoLegend() + 
  #FontSize(x.title = 13, y.title = 13)
  theme(axis.title = element_text(size=13)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(axis.title.y = element_text(size = 13, vjust = -3), 
        axis.text.y = element_text(size = 10, margin = margin(0,1,0,0))) +
  theme(axis.title.x = element_text(size = 13, vjust = 1), 
        axis.text.x = element_text(size = 10, margin = margin(1,0,1,0))) +
  theme(strip.text = element_blank()) +
  scale_x_continuous(breaks=seq(-10,10,10)) +
  scale_shape_manual(values=16)
dev.off()


#======================================================================================
# Save
#======================================================================================
saveRDS(sub.o, paste0("data_09032020/svz_subset_OSKM_",Sys.Date(),".rds"))


#=======================================================================================
# Differential expression by MAST for each cell type
#=======================================================================================

sub.o <- readRDS(paste0("data_09032020/svz_subset_OSKM_2020-09-24.rds"))
sub.oo <- subset(sub.o, age=="old")

## Append treatment to celltype label (low res) and add to metadata
cellAge <- paste0(sub.oo@meta.data$Celltype, "_", sub.oo@meta.data$Treatment)
m <- data.frame("Celltype_Treatment" = cellAge)
rownames(m) <- rownames(sub.oo@meta.data)
sub.oo <- AddMetaData(object = sub.oo, metadata = m)

# Assign as main identity
Idents(sub.oo) <- sub.oo@meta.data$Celltype_Treatment

# Assign default assay
DefaultAssay(sub.oo) <- "RNA"

data.frame(unclass(table(sub.o$Celltype, sub.oo$Treatment)))


#-----------------------------------
# Old+OSKM vs. Old control
#-----------------------------------

# Exclude celltypes with less than 3 cells per condition and unknown clusters
CELLTYPES <- c("Astrocyte_qNSC", "Oligodendrocyte", "Microglia", "Neuroblast", "aNSC_NPC", "Endothelial", "Mural", "Macrophage")
mast_list <- vector(mode="list", length = 8)
names(mast_list) = CELLTYPES

# Compare 2Dox0 and untr for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  svz.de.2Dox0.untr <- FindMarkers(object = sub.oo,
                             ident.1 = paste0(CELLTYPE, "_2Dox0"),
                             ident.2 = paste0(CELLTYPE, "_untr"),
                             assay = "RNA",
                             slot = "data",
                             only.pos = FALSE, 
                             min.pct = 0,
                             logfc.threshold = 0,
                             test.use = "MAST")
  svz.de.2Dox0.untr <- as.data.frame(svz.de.2Dox0.untr)
  svz.de.2Dox0.untr <- rownames_to_column(svz.de.2Dox0.untr, var = "gene")
  svz.de.2Dox0.untr$celltype <- rep(CELLTYPE, dim(svz.de.2Dox0.untr)[1])
  mast_list[[CELLTYPE]] <- svz.de.2Dox0.untr
}

# Combine all matrices into one dataframe
mast_df_2Dox0.untr <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df_2Dox0.untr  <- rbind(mast_df_2Dox0.untr, mast_list[[CELLTYPE]])
}
dim(mast_df_2Dox0.untr) # 133816 by 7
head(mast_df_2Dox0.untr) 

# Calculate adjusted p-values using FDR, and add to dataframe                            
mast_df_2Dox0.untr$p_adj_fdr <- p.adjust(mast_df_2Dox0.untr$p_val, method = "fdr")
head(mast_df_2Dox0.untr)

# Save differential expression results
write.csv(mast_df_2Dox0.untr, "data_MAST/de.mast.OSKM-old.2Dox0.untr_df.csv")

for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  write.csv(mast_list[[CELLTYPE]], paste0("data_MAST/de.mast.OSKM-old.2Dox0.untr_df_",CELLTYPE,".csv"))
}


#-----------------------------------
# Old control vs. young control
#-----------------------------------

sub.ou <- subset(sub.o, Treatment == "untr")
DefaultAssay(sub.ou) <- "RNA"

# Assign as main identity
Idents(sub.ou) <- sub.ou@meta.data$Celltype_Age

#Check cell counts
data.frame(unclass(table(sub.ou$Celltype, sub.ou$age)))

CELLTYPES <- unique(sub.ou@meta.data$Celltype.LowRes)
CELLTYPES <- c("Astrocyte_qNSC", "Oligodendrocyte", "OPC","Microglia", "Neuroblast", "aNSC_NPC", "Endothelial", "Ependymal","Mural", "Macrophage")
mast_list <- vector(mode="list", length = 10)
names(mast_list) = CELLTYPES

# Compare young and old for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  svz.de.old.young <- FindMarkers(object = sub.ou,
                                   ident.1 = paste0(CELLTYPE, "_old"),
                                   ident.2 = paste0(CELLTYPE, "_young"),
                                   only.pos = FALSE, 
                                   min.pct = 0,
                                   logfc.threshold = 0,
                                   test.use = "MAST")
  svz.de.old.young <- as.data.frame(svz.de.old.young)
  svz.de.old.young <- rownames_to_column(svz.de.old.young, var = "gene")
  svz.de.old.young$celltype <- rep(CELLTYPE, dim(svz.de.old.young)[1])
  mast_list[[CELLTYPE]] <- svz.de.old.young
}

# Combine all matrices into one dataframe
mast_df_old.young <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df_old.young  <- rbind(mast_df_old.young, mast_list[[CELLTYPE]])
}
dim(mast_df_old.young) # 166740 by 7
head(mast_df_old.young) 

# Calculate adjusted p-values using FDR, and add to dataframe                            
mast_df_old.young$p_adj_fdr <- p.adjust(mast_df_old.young$p_val, method = "fdr")
head(mast_df_old.young)

# Save differential expression results
write.csv(mast_df_old.young, "data_09032020/de.mast.OSKM-untr.old.young_df.csv")



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