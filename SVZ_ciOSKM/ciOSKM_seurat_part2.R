# SVZ-targeted partial reprogramming (c+iOSKM))
# LMO + SVZ Seurat Processing
# 2 cohorts across 6 total 10x lanes

# PART TWO

# v3.2.2 Seurat

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

setwd("~/Dropbox/10x_expAC/3.Seurat") # LOCAL

#======================================================================================
# Load data
markers <- readRDS("data/svz.combined.sing.MULTI.CLR.2.rescued.markers.2022-05-20.rds")
svz_sing <- readRDS("data/svz.combined.sing.MULTI.CLR.2.rescued.2022-05-20.rds")
t20 <- tableau_color_pal(palette = "Tableau 20")(20)
t10 <- tableau_color_pal(palette = "Tableau 10")(10)
t4 <- tableau_color_pal(palette = "Tableau 10")(4)[c(2,3,4,1)]


svz_sing <- RunUMAP(svz_sing, dims = 1:14, min.dist = .3, spread = .75, seed.use = 22)
DimPlot(svz_sing, group.by = "seurat_clusters", label = TRUE, label.size = 5, pt.size = .6)

ggsave(paste0("plots2/umap.seed", 22, ".pdf"), width=7, height=6)


#======================================================================================
#===== Cell type markers =====#
#======================================================================================

# Microglia, olig, and endothelial marker genes
genes <- c("C1qc", "Junb", "Cd52", "Opalin", "Tspan2", "Klk6", "Cxcl12", "Ly6c1", "Itm2a")
FeaturePlot(svz_sing, features = genes)
ggsave("plots2/umap.micro-olig-endo-markers.pdf", width = 12, height = 10)


# NSC lineage 
genes <- c("Clu", "Aldoc", "Pclaf", "Cenpf","Meg3", "Tubb3")
FeaturePlot(svz_sing, features = genes)
ggsave("plots2/umap.NSClineage-markers.pdf", width = 8, height = 10)


# T cell markers
# Cd3, Cd8, Cd4, interferon, activation, and tissue retention genes
cd.genes <- grep(pattern = "^Cd3", x = rownames(x = svz_sing$SCT@data), value = TRUE)
genes <- c("Cd3g", "Cd3d", "Cd3e", "Cd8a", "Cd4", "Ifng", "Cd69", "Xcl1", "Itga4") #"Cd8b"
FeaturePlot(svz_sing, features = genes, order = T)
ggsave("plots2/umap.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)

# Check full dataset including doublets and negatives
FeaturePlot(svz, features = genes, order = F)
ggsave("plots2/umap.doub.neg.sing.Tcell-Cd3-Cd8-Cd4-Ifng-Cd69-Xcl1-Itga4.pdf", width = 12, height = 10)
genes <- c("Cd3e", "Cd8a")
FeaturePlot(svz, features = genes, split.by = "MULTI_classification.global", order = TRUE, pt.size = 0.8)
ggsave("plots2/umap.split.doub.neg.sing.Cd3-Cd8.pdf", width = 12, height = 8)

FeaturePlot(svz_sing, features = genes, split.by = "orig.ident", order = TRUE, pt.size = 0.8)
ggsave("plots2/umap.splitlane.doub.neg.sing.Cd3-Cd8.pdf", width = 20, height = 8)


# Mural: vascular smooth muscle cells (Acta2, Myl9, Pdlim3) and pericytes (Kcnj8, Ifitm1,Rgs5)
# Other marker genes from https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
FeaturePlot(svz_sing, features = c("Pdlim3","Acta2","Myl9","Kcnj8","Rgs5" ,"Ifitm1","Pdgfrb"),order=T)
ggsave("plots2/umap.mural.pdf", height = 9, width = 10)

# OPCs. some Pdgfra+. Opalin-Mog-.
FeaturePlot(svz_sing, features = c("Pdgfra","Opalin","Mbp","Mog"),order=T)
ggsave("plots2/umap.OPC.pdf", height = 7, width = 8)


# Microglia/macrophage subtypes
# Cx3cr1 is a marker for interstitial microglia, Mrc1 is a marker for perivascular microglia and macrophages. https://www-nature-com.stanford.idm.oclc.org/articles/nature25739 
FeaturePlot(svz_sing, features = c("Cx3cr1", "Mrc1","Ifi27l2a", "Lyz2"),order=T)
ggsave("plots2/umap.microglia.pdf", height = 6, width = 6)


# Monocyte/macrophage-specific markers from https://doi-org.stanford.idm.oclc.org/10.1016/j.immuni.2018.11.004 
FeaturePlot(svz_sing,features=c("F13a1","H2-Aa","Ccr2","Lyve1","Mgl2", "P2ry12"), order=T)
ggsave("plots2/umap.macrophage-Stevens.pdf", height = 9, width = 6)


# Ependymal cells
# aSMA (acta2) is used as a marker to sort. Other ependymal markers are also specific to this cluster. https://doi-org.stanford.idm.oclc.org/10.1016/j.cell.2018.03.063 
FeaturePlot(svz_sing, features = c("Acta2", "Rarres2", "Tmem212", "Mia"), order=T)
ggsave("plots2/umap.ependymal.pdf", height = 7, width = 8)


# SVZ lineage markers
FeaturePlot(svz_sing, features = c("Gfap", "Slc1a3", "Egfr", "Nes",  "Prom1", "Ncam1"), order = T)
ggsave("plots2/umap.NSCmarkers.Gfap.Slc1a3.Egfr.Nes.Prom1.Ncam1.pdf", width = 8, height = 10)


# (putative neurons)
FeaturePlot(svz_sing, features = c("Neurod1", "Neurod4", "Snap25","Chgb"), order=T) & NoLegend() & NoAxes()
ggsave("plots2/umap.neuron.nolegend.noaxes.pdf", width = 6, height = 6)


# (putative OPCs)
FeaturePlot(svz_sing, features = c("Bmp4", "Neu4", "Tnr","C1ql1"), order = T)



#======================================================================================
# Assign celltypes
#======================================================================================
#25 clusters with resolution = 0.5
new.cluster.ids <- c("Microglia_1",
                     "Microglia_2",
                     "Oligodendrocyte_1",
                     "Endothelial_1",
                     "Neuroblast_1",
                     "Mono_Mac_1",
                     "Astrocyte_qNSC_1",
                     "Oligodendrocyte_2",
                     "aNSC_NPC_1",
                     "Neuroblast_2",
                     "Microglia_3",
                     "T-cell",
                     "Endothelial_2",
                     "Mural.Pericyte",
                     "Mono_Mac_2",
                     "Mono_Mac_3",
                     "Ependymal",
                     "Mono_Mac_4",
                     "Mural.VSMC",
                     "Neuron",
                     "Mono_Mac_5",
                     "Astrocyte_qNSC_2",
                     "OPC",
                     "Oligodendrocyte_3",
                     "aNSC_NPC_2")

Idents(svz_sing) <- "seurat_clusters"
names(new.cluster.ids) <- levels(svz_sing)
svz_sing <- RenameIdents(svz_sing, new.cluster.ids)
svz_sing[["Celltype.25"]] <- Idents(svz_sing)
unique(svz_sing@meta.data$Celltype.25)


#======================================
# Refine cell type annotations

# Small qNSC tail that is actually OPCs based on Pdgfra expression ----
plot <- DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.25", label = TRUE, label.size = 5, pt.size = .6)
select.cells <- CellSelector(plot = plot) #74 cells
head(select.cells)
Idents(svz_sing) <- "Celltype.25"
Idents(svz_sing, cells = select.cells) <- "sub"
newcells.markers <- FindMarkers(svz_sing, ident.1 = "sub", ident.2 = "Astrocyte_qNSC_1", min.diff.pct = 0.3,
                                only.pos = FALSE)
head(newcells.markers)
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# Pdgfra 3.062867e-162   2.714330 0.875 0.005 5.911026e-158
# Lhfpl3 5.716528e-158   2.817229 0.891 0.008 1.103233e-153
# Opcml  1.452829e-157   2.188874 0.875 0.006 2.803814e-153
# Tnr    1.534151e-153   2.507127 0.875 0.008 2.960758e-149
# Gjc3   2.673267e-148   2.326924 0.875 0.011 5.159138e-144
# Cspg4  4.849796e-135   1.979222 0.797 0.009 9.359622e-131

# Correct annotation
Idents(svz_sing, cells = select.cells) <- "OPC"
svz_sing[["Celltype.25"]] <- Idents(svz_sing)
DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.25", label = TRUE, label.size = 5, pt.size = .6)


# Small cluster of putative aNSC/NPC that is separate in UMAP space ----
plot <- DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.25", label = TRUE, label.size = 5, pt.size = .6)
select.cells <- CellSelector(plot = plot) #47 cells
head(select.cells)
Idents(svz_sing) <- "Celltype.25"
Idents(svz_sing, cells = select.cells) <- "sub"
newcells.markers <- FindMarkers(svz_sing, ident.1 = "sub", ident.2 = "aNSC_NPC_1", min.diff.pct = 0.3,
                                only.pos = FALSE)
head(newcells.markers)
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# Ptgs1  5.690952e-139   1.892813 0.936 0.013 1.098297e-134
# Pld4   4.845828e-133   1.818276 0.957 0.018 9.351964e-129
# Anxa3  1.460701e-132   2.027077 0.915 0.014 2.819007e-128
# Ly86   1.224406e-130   4.322215 1.000 0.026 2.362981e-126
# Trem2  1.193408e-128   3.672497 0.979 0.024 2.303159e-124
# Olfml3 5.399357e-128   2.617804 0.936 0.019 1.042022e-123

# Could be doublets with microglia based on expression of microglia markers? Strong proliferation signature. 


#======================================

Idents(svz_sing) <- "Celltype.25"
svz_sing@meta.data$Celltype.25 <- Idents(svz_sing)
levels(unique(svz_sing@meta.data$Celltype.25))
table(svz_sing@meta.data$Celltype.25)

# Assign lower resolution celltypes
# order changed (OPC inserted first)
low.res.ids <-  c("OPC",
                  "Microglia",
                  "Microglia",
                  "Oligodendrocyte",
                  "Endothelial",
                  "Neuroblast",
                  "Mono_Mac",
                  "Astrocyte_qNSC",
                  "Oligodendrocyte",
                  "aNSC_NPC",
                  "Neuroblast",
                  "Microglia",
                  "T-cell",
                  "Endothelial",
                  "Pericyte",
                  "Mono_Mac",
                  "Mono_Mac",
                  "Ependymal",
                  "Mono_Mac",
                  "VascSmoothMuscle",
                  "Neuron",
                  "Mono_Mac",
                  "Astrocyte_qNSC",
                  "Oligodendrocyte",
                  "aNSC_NPC")

Idents(svz_sing) <- svz_sing@meta.data$Celltype.25
names(low.res.ids) <- levels(svz_sing)
svz_sing <- RenameIdents(svz_sing, low.res.ids)
svz_sing[["Celltype.LowRes"]] <- Idents(svz_sing)
unique(svz_sing@meta.data$Celltype.LowRes)
Idents(svz_sing) <- svz_sing@meta.data$Celltype.LowRes


DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.25", label = TRUE, label.size = 4, repel = T, pt.size = .6)
ggsave("plots2/umap.celltype.label.pdf", width=9, height=6)

DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.LowRes", label = TRUE, label.size = 4, repel = T, pt.size = .6, cols = t20)
ggsave("plots2/umap.celltype.lowres.label.pdf", width=7.5, height=6)

DimPlot(svz_sing, reduction = "umap", group.by = "Celltype.LowRes", label = F, pt.size = .6, cols = t20)
ggsave("plots2/umap.celltype.lowres.pdf", width=7.5, height=6)



# DEFINE COLORS BY CELLTYPE --------------------
CELLSALL <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuron",
           "Oligodendrocyte", "OPC", "Endothelial",
           "Microglia", "Mono_Mac", "T-cell", 
           "Pericyte", "VascSmoothMuscle", "Ependymal")
svz_sing$Celltype <- factor(svz_sing$Celltype.LowRes, levels = CELLSALL, ordered = T)

cellcolors <- c("#03c03c",  "#4E79A7", "#966fd6", "#499894",
                "#ffdf00", "#B6992D", "#FFBE7D", 
                "#aec6cf", "#86BCB6", "#f78d76", 
                "#db7093", "#A52A2A", "#2F4F4F")
cellcolors <- alpha(cellcolors, 0.6)
names(cellcolors) <- CELLSALL
names(cellcolors) <- levels(svz_sing$Celltype)
# ----------------------------------------


DimPlot(svz_sing, group.by = "Celltype", label = TRUE, cols = cellcolors, label.size = 4, pt.size = .6, repel = T) + labs(title=NULL) 
ggsave(paste0("plots2/umap.celltypes.labels.pdf"), width=7, height=6)

DimPlot(svz_sing, group.by = "Celltype", cols = cellcolors, pt.size = 1, shape.by = "orig.ident") + labs(title=NULL) + scale_shape_manual(values=c(16,16,16,16,16,16), guide="none")
ggsave(paste0("plots2/umap.celltypes.pdf"), width=8, height=6)



#====================================
# Add sample information
#====================================
unique(svz_sing$MULTI_classification_rescued)
age_treatment <- plyr::mapvalues(x = svz_sing$MULTI_classification_rescued, 
                             from = c("A-1108L-3.9-M-lane1-untr-TGTGATGG", 
                                      "B-1110L-3.9-M-lane1-untr-TCAATGGC", 
                                      "C-L0111-28.1-M-lane2-untr-CTCTAGAC", 
                                      "D-L0112-28.1-M-lane3-Dox-ACCAATGC", 
                                      "E-L0114-28.1-M-lane2-untr-AGTTGCGT", 
                                      "F-L0116-28.1-F-lane3-Dox-CGAACAAG",
                                      "G-L0117-28.1-F-lane2-untr-GTACCTGT",
                                      "H-L125-27.4-M-lane3-Dox-GAAGCTTG",
                                      "I-1114L-4.0-F-lane4-untr-AAGTACGC",
                                      "J-1115L-4.0-F-lane4-untr-ATTCGCAC",
                                      "K-1116L-4.0-F-lane4-untr-GAGTCGAT",
                                      "L-L131-27.4-F-lane6-Dox-AAGGCTAG",
                                      "M-L126-27.4-M-lane5-untr-CAGTTAGG",
                                      "N-L128-27.4-F-lane6-Dox-AACCGAAC",
                                      "O-L132-27.4-F-lane5-untr-AAGCAGTC",
                                      "P-L133-26.3-M-lane6-Dox-GAATCAGG",
                                      "Q-L0108-28.8-M-lane5-untr-ACTCGAAG"), 
                             to = c("young_untr", 
                                    "young_untr", 
                                    "old_untr", 
                                    "old_Dox", 
                                    "old_untr", 
                                    "old_Dox", 
                                    "old_untr", 
                                    "old_Dox", 
                                    "young_untr", 
                                    "young_untr",
                                    "young_untr", 
                                    "old_Dox", 
                                    "old_untr",
                                    "old_Dox", 
                                    "old_untr",
                                    "old_Dox", 
                                    "old_untr"))
svz_sing$age_treatment <- age_treatment


age <- plyr::mapvalues(x = svz_sing$age_treatment, 
                                from = c("young_untr", 
                                         "old_untr", 
                                         "old_Dox"),
                                to = c("young", 
                                       "old", 
                                       "old"))
svz_sing$age <- age

treatment <- plyr::mapvalues(x = svz_sing$age_treatment, 
                       from = c("young_untr", 
                                "old_untr", 
                                "old_Dox"),
                       to = c("untr", 
                              "untr", 
                              "Dox"))

svz_sing$treatment <- treatment


# Adjust factors-----------------------------
AGE <- c("young", "old")
svz_sing$age <- factor(svz_sing$age, levels=AGE, ordered=T)
table(svz_sing[[c("MULTI_classification_rescued", "age")]])

TREATMENT <- c("untr", "Dox")
svz_sing$treatment <- factor(svz_sing$treatment, levels=TREATMENT, ordered=T)
table(svz_sing[[c("MULTI_classification_rescued", "treatment")]])

AGE.TREATMENT <- c("young_untr", "old_untr", "old_Dox")
svz_sing$age_treatment <- factor(svz_sing$age_treatment, levels=AGE.TREATMENT, ordered=T)
table(svz_sing[[c("MULTI_classification_rescued", "age_treatment")]])


# SET COLORS FOR TREATMENTS-----------------------------
doxcolors <- c("turquoise","tomato", "#612a95")
names(doxcolors) <- c("young_untr","old_untr","old_Dox")
# ---------------------------------------------------------

Idents(svz_sing) <- "age_treatment"
DimPlot(svz_sing, reduction = "umap", label = F, label.size = 5, pt.size = 1, shape.by="age", cols = alpha(doxcolors,0.4)) + scale_shape_manual(values=c(16,16), guide = "none")
ggsave("plots2/umap.age.treatment.pdf", width=7.5, height=6)


#======================================================================================
# SAVE ----
#======================================================================================

saveRDS(svz_sing, paste0("data/svz_celltypes_metadata_", Sys.Date(), ".rds"))
svz_sing <- readRDS(paste0("data/svz_celltypes_metadata_", "2022-05-20", ".rds"))

#======================================================================================

# Marker genes after annotation for table
Idents(svz_sing) <- svz_sing@meta.data$Celltype
newmarkers <- FindAllMarkers(object=svz_sing)
write.csv(newmarkers, paste0('data/svz.combined.sing.MULTI.CLR.2.rescued.markers.annotated', Sys.Date(), '.csv'))

# Table of number of each cell type in each sample
num <- data.frame(unclass(table(svz_sing$Celltype, svz_sing$MULTI_classification_rescued)))
write.csv(num,paste0("data/table.celltype-number.MULTI-ID-rescued.",Sys.Date(),".csv"))



#======================================================================================
# Downsample by condition
#======================================================================================

table(svz_sing$age_treatment)
#Downsample to lowest number of cells per treatment condition
# young_untr   old_untr    old_Dox 
# 9915       3112       3264 

Idents(svz_sing) <- 'age_treatment'
downsampled <- subset(svz_sing, downsample = 3112)

table(downsampled$age_treatment)
# young_untr   old_untr    old_Dox 
# 3112       3112       3112 

downsampled$age_treatment <- factor(downsampled$age_treatment, levels=c("young_untr","old_untr","old_Dox"), ordered=T)

DimPlot(downsampled, reduction = "umap", label = F, label.size = 5, pt.size = .8, split.by = "age_treatment", cols = alpha(doxcolors,0.5), shape.by = "age") + NoLegend() +scale_shape_manual(values=c(16,16), guide = "none")
ggsave("plots2/umap.downsampled.split-age-treatment.nolegend.pdf", width=11, height=4.5)



#=======================================================================================
# Differential Expression for each cell type
#=======================================================================================

setwd("~/Dropbox/10x_expAC/3.Seurat")
svz_sing <- readRDS("data/svz_celltypes_metadata_2022-05-20.rds")

## Append treatment to celltype label (low res) and add to metadata
cellAge <- paste0(svz_sing@meta.data$Celltype, "_", svz_sing@meta.data$age_treatment)
m <- data.frame("Celltype_Treatment" = cellAge)
rownames(m) <- rownames(svz_sing@meta.data)
svz_sing <- AddMetaData(object = svz_sing, metadata = m)

# Assign as main identity
Idents(svz_sing) <- svz_sing@meta.data$Celltype_Treatment

# Assign default assay
DefaultAssay(svz_sing) <- "RNA"


data.frame(unclass(table(svz_sing$Celltype, svz_sing$age_treatment)))
#                  young_untr old_untr old_Dox
# Astrocyte_qNSC          654      199      81
# aNSC_NPC                591       86     149
# Neuroblast             1338      240     297
# Neuron                    0        0     194
# Oligodendrocyte        1154      679     398
# OPC                      69       35      20
# Endothelial             945      576     352
# Microglia              2586      986     978
# Mono_Mac               1722       64     357
# T-cell                  501       47      94
# Pericyte                238      129     101
# VascSmoothMuscle        104       49      46
# Ependymal                13       22     197


CELLTYPES <- unique(svz_sing@meta.data$Celltype.LowRes)
mast_list <- vector(mode="list", length = 13)
names(mast_list) = CELLTYPES


#====
# old untr vs. young untr
#====

# Exclude celltypes with less than 3 cells per condition (exclude neurons- only found in old+Dox condition)
CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Endothelial", "Oligodendrocyte", "OPC", "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
mast_list <- vector(mode="list", length = 12)
names(mast_list) = CELLTYPES

# Compare old and young untr for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  svz.de.old.young <- FindMarkers(object = svz_sing,
                             ident.1 = paste0(CELLTYPE, "_old_untr"),
                             ident.2 = paste0(CELLTYPE, "_young_untr"),
                             assay = "RNA",
                             slot = "data",
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
dim(mast_df_old.young) # 372708 by 7
head(mast_df_old.young) 


# Calculate adjusted p-values using FDR, and add to dataframe                            
mast_df_old.young$p_adj_fdr <- p.adjust(mast_df_old.young$p_val, method = "fdr")
head(mast_df_old.young)

# Save differential expression results
write.csv(mast_df_old.young, "data_MAST/de.mast_old.young_df.csv")

for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  write.csv(mast_list[[CELLTYPE]], paste0("data_MAST/de.mast_old.young_df_",CELLTYPE,".csv"))
}


#====
# old Dox vs. old untr
#====

# Exclude celltypes with less than 3 cells per condition (exclude neurons- only found in old+Dox condition)
CELLTYPES <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Endothelial", "Oligodendrocyte", "OPC", "Microglia", "Mono_Mac", "T-cell", "Pericyte", "VascSmoothMuscle", "Ependymal")
mast_list <- vector(mode="list", length = 12)
names(mast_list) = CELLTYPES


# Compare old Dox vs. old untr for each celltype, save each matrix in a list. 
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  svz.de.Dox.untr <- FindMarkers(object = svz_sing,
                                   ident.1 = paste0(CELLTYPE, "_old_Dox"),
                                   ident.2 = paste0(CELLTYPE, "_old_untr"),
                                   assay = "RNA",
                                   slot = "data",
                                  only.pos = FALSE,
                                  min.pct = 0,
                                  logfc.threshold = 0,
                                  test.use = "MAST")
  
  svz.de.Dox.untr <- as.data.frame(svz.de.Dox.untr)
  svz.de.Dox.untr <- rownames_to_column(svz.de.Dox.untr, var = "gene")
  svz.de.Dox.untr$celltype <- rep(CELLTYPE, dim(svz.de.Dox.untr)[1])
  mast_list[[CELLTYPE]] <- svz.de.Dox.untr
}

# Combine all matrices into one dataframe
mast_df_Dox.untr <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df_Dox.untr  <- rbind(mast_df_Dox.untr, mast_list[[CELLTYPE]])
}
dim(mast_df_Dox.untr) # 372708 by 7
head(mast_df_Dox.untr) 


# Calculate adjusted p-values using FDR, and add to dataframe                            
mast_df_Dox.untr$p_adj_fdr <- p.adjust(mast_df_Dox.untr$p_val, method = "fdr")
head(mast_df_Dox.untr)


# Save differential expression results
write.csv(mast_df_Dox.untr, "data_MAST/de.mast_oldDox.olduntr_df.csv")

for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  write.csv(mast_list[[CELLTYPE]], paste0("data_MAST/de.mast_oldDox.olduntr_df_",CELLTYPE,".csv"))
}


