# Run GSEA analysis on gene list from DESeq2 analysis

library(tidyverse)
library(biomaRt)
library(fgsea)
library(DT)
library(org.Mm.eg.db)
library(msigdbr)
library(NCmisc)
library(clusterProfiler)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())


setwd("~/Dropbox/bulkRNAseq_expAL/2.GSEA/")

# Load pathways. Genes are as human gene symbols.
pathways.go <- gmtPathways("~/Dropbox/10x_OSKM_2/4.GSEA/pathways/c5.all.v7.1.symbols.gmt") #GO gene sets


# ====================================
# A18h old Dox vs. untr ----
# ====================================

COMPARISON <- "A18h-old-Dox_A18h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_A18h-old-Dox_A18h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>% 
  mutate_all(na_if,"") %>%
  na.omit()
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
    na.omit() %>%
    distinct() %>%
    group_by(hsapiens_homolog_associated_gene_name) %>% 
    dplyr::summarize(stat=mean(z)) %>% 
    arrange(desc(stat))
  
  # Genes ranked by Z score (not log fold change)
  ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))
  
fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
print(head(fgseaResTidy))
  
fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))
  
fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# Plotting ----

data <- fgseaResTidy 

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^...", "", data$pathway)

ggplot(data %>% dplyr::filter(padj<0.05) %>% dplyr::filter(dense_rank(NES)<=15 | dense_rank(-NES)<=15), aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots/",COMPARISON,"_fgsea_GO_dotplot_padj0.05_NEStop15bottom15.pdf"), width=6, height=4)  


# ====================================
# A18h old vs. young ----
# ====================================

COMPARISON <- "A18h-old-untr_A18h-young-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_A18h-old-untr_A18h-young-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D06h old Dox vs. untr ----
# ====================================

COMPARISON <- "D06h-old-Dox_D06h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D06h-old-Dox_D06h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D06h old vs. young ----
# ====================================

COMPARISON <- "D06h-old-untr_D06h-young-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D06h-old-untr_D06h-young-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D18h old Dox vs. untr ----
# ====================================

COMPARISON <- "D18h-old-Dox_D18h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D18h-old-Dox_D18h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D18h old vs. young ----
# ====================================

COMPARISON <- "D18h-old-untr_D18h-young-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D18h-old-untr_D18h-young-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D48h old Dox vs. untr ----
# ====================================

COMPARISON <- "D48h-old-Dox_D48h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D48h-old-Dox_D48h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# D48h old vs. young ----
# ====================================

COMPARISON <- "D48h-old-untr_D48h-young-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_D48h-old-untr_D48h-young-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")



# ====================================
# A06h old Dox vs. A18h old untr ----
# ====================================

COMPARISON <- "A06h-old-Dox_A18h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_A06h-old-Dox_A18h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")


# Plotting ----

data <- fgseaResTidy 

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^...", "", data$pathway)

ggplot(data %>% dplyr::filter(padj<0.05) %>% dplyr::filter(dense_rank(NES)<=15 | dense_rank(-NES)<=15), aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots/",COMPARISON,"_fgsea_GO_dotplot_padj0.05_NEStop15bottom15.pdf"), width=6.5, height=4)  


# ====================================
# A48h old Dox vs. A18h old untr ----
# ====================================

COMPARISON <- "A48h-old-Dox_A18h-old-untr"

# Load cell type differential expression results
mast <- read.csv("../1.DESeq2/data/DESeqRes_designCondition_A48h-old-Dox_A18h-old-untr.csv")
mast$z <- p.to.Z(mast$pvalue) * sign(mast$log2FoldChange)
mast$z.adj <- p.to.Z(mast$padj) * sign(mast$log2FoldChange)

# Get extra gene information, combine with DE results
df <- inner_join(mast, bm, by = c("Gene" = "ensembl_gene_id"))

# Rank list
df <- df %>% dplyr::select(hsapiens_homolog_associated_gene_name, z) %>% 
  na.omit() %>%
  distinct() %>%
  group_by(hsapiens_homolog_associated_gene_name) %>% 
  dplyr::summarize(stat=mean(z)) %>% 
  arrange(desc(stat))

# Genes ranked by Z score (not log fold change)
ranks <- deframe(df)

# Run the Gene Set Enrichment Analysis
fgseaRes <- fgsea(pathways=pathways.go, stats=ranks) 
print(head(fgseaRes))

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
print(head(fgseaResTidy))

fgseaResTidy <- as.data.frame(fgseaResTidy)

fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste, collapse = ", ", character(1L))

fname <- paste0("data/", COMPARISON, "_fgsea_GO_leadingedge_", Sys.Date(),".txt")
write.table(fgseaResTidy, file=fname, sep = "\t")


# Plotting ----

data <- fgseaResTidy 

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^..", "", data$pathway)

ggplot(data %>% dplyr::filter(padj<0.05) %>% dplyr::filter(dense_rank(NES)<=15 | dense_rank(-NES)<=15), aes(x=NES, y=reorder(pathway,NES), size=size, color=padj)) + 
  geom_point(alpha = 1) + 
  scale_color_viridis(limits=c(0,0.055),oob = scales::squish)+ #, guide='none') +
  theme_minimal() +
  labs(x="Normalized Enrichment Score", #y="Pathway",title="GO terms from GSEA",
       subtitle=paste0(COMPARISON)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=10), 
        legend.text = element_text(size=10)) 

ggsave(paste0("plots/",COMPARISON,"_fgsea_GO_dotplot_padj0.05_NEStop15bottom15.pdf"), width=5.5, height=4)  



# =======================
# Scatter plots ----
# =======================

# Load data and combine ----
# A18h
gsea_yo <- read.table("data/A18h-old-untr_A18h-young-untr_fgsea_GO_leadingedge_2023-05-31.txt")
gsea_dox <-read.table("data/A18h-old-Dox_A18h-old-untr_fgsea_GO_leadingedge_2023-05-30.txt")

gsea_yo$timepoint <- "A18h"
gsea_dox$timepoint <- "A18h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_yo, gsea_dox)

# D06h
gsea_yo <- read.table("data/D06h-old-untr_D06h-young-untr_fgsea_GO_leadingedge_2023-06-01.txt")
gsea_dox <-read.table("data/D06h-old-Dox_D06h-old-untr_fgsea_GO_leadingedge_2023-06-01.txt")

gsea_yo$timepoint <- "D06h"
gsea_dox$timepoint <- "D06h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_complete, gsea_yo, gsea_dox)

# D18h
gsea_yo <- read.table("data/D18h-old-untr_D18h-young-untr_fgsea_GO_leadingedge_2023-06-01.txt")
gsea_dox <-read.table("data/D18h-old-Dox_D18h-old-untr_fgsea_GO_leadingedge_2023-06-01.txt")

gsea_yo$timepoint <- "D18h"
gsea_dox$timepoint <- "D18h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_complete, gsea_yo, gsea_dox)

# D48h
gsea_yo <- read.table("data/D48h-old-untr_D48h-young-untr_fgsea_GO_leadingedge_2023-06-01.txt")
gsea_dox <-read.table("data/D48h-old-Dox_D48h-old-untr_fgsea_GO_leadingedge_2023-06-01.txt")

gsea_yo$timepoint <- "D48h"
gsea_dox$timepoint <- "D48h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete <- bind_rows(gsea_complete, gsea_yo, gsea_dox)

# A06h
gsea_yo <- read.table("data/A18h-old-untr_A18h-young-untr_fgsea_GO_leadingedge_2023-05-31.txt") #young vs old comparisons are from A18h timepoint
gsea_dox <-read.table("data/A06h-old-Dox_A18h-old-untr_fgsea_GO_leadingedge_2023-06-02.txt") #old untr reference is A18h

gsea_yo$timepoint <- "A06h"
gsea_dox$timepoint <- "A06h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete2 <- bind_rows(gsea_complete, gsea_yo, gsea_dox)

# A48h
gsea_yo <- read.table("data/A18h-old-untr_A18h-young-untr_fgsea_GO_leadingedge_2023-05-31.txt") #young vs old comparisons are from A18h timepoint
gsea_dox <-read.table("data/A48h-old-Dox_A18h-old-untr_fgsea_GO_leadingedge_2023-06-02.txt") #old untr reference is A18h

gsea_yo$timepoint <- "A48h"
gsea_dox$timepoint <- "A48h"

gsea_yo$comparison <- "yo"
gsea_dox$comparison <- "dox"

gsea_complete2 <- bind_rows(gsea_complete2, gsea_yo, gsea_dox)


# Save ---
dim(gsea_complete2) #122220    10
gsea_complete2$comparison <- factor(gsea_complete2$comparison, levels=c("yo","dox"), ordered=T)
write.table(gsea_complete2, paste0("data/GSEA_complete_allAct_", Sys.Date(), ".txt"))

# Load ----
gsea_complete2 <- read.table("data/GSEA_complete_allAct_2023-06-02.txt")



# Restructure for plotting ----

widegsea2 <- reshape(gsea_complete2, idvar = c("pathway","timepoint"), timevar = "comparison", v.names=c("NES","padj","pval","log2err","size"), direction="wide")
widegsea2 <- subset(widegsea2, select=-c(ES,leadingEdge))


# Set significance factors ----
widegsea2$a <- "1"
a <- widegsea2$a
a[widegsea2$padj.dox<0.05 & widegsea2$padj.yo<0.05] <- "2" #Significant in both
widegsea2$a <- a

sigalpha <- c(0.2,0.9)
names(sigalpha) <- c("1","2")

sigshape <- c(1,19)
names(sigshape) <- c("1","2")

timecols <- c("#000004FF", "#420A68FF", "#932667FF", "#DD513AFF", "#FCA50AFF", "#FDE725FF")
names(timecols) <- c("A06h", "A18h", "A48h", "D06h", "D18h", "D48h")


# Plot ----
ggplot(widegsea2 %>% filter(padj.yo<0.05|padj.dox<0.05), aes(x=NES.yo, y=NES.dox, color=timepoint, alpha=a)) + 
  geom_point(size=0.8, aes(alpha=a, shape=a)) + 
  facet_wrap(~timepoint, nrow=2) +
  scale_color_manual(values = timecols, guide='none') +
  scale_alpha_manual(values = sigalpha, guide='none') +
  scale_shape_manual(values = sigshape, guide='none') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 9), 
        axis.text.x = element_text(size=9), axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        panel.border = element_rect(color = "darkgray", size=1, fill=NA)) +
  geom_vline(xintercept=0) + geom_hline(yintercept=0) 

ggsave("plots/scatterGSEA_facet-all_Dox-untr_v_old-young_subset-padj0.05-either_fill-padj0.05-both.pdf", height = 4.5, width = 6)




#==========
# Dotplots- Compare Dox-v-old and old-v-young, across timepoints ----
#==========

setwd("~/Dropbox/bulkRNAseq_expAL/2.GSEA/")

# Load data and combine
gsea_complete2 <- read.table("data/GSEA_complete_allAct_2023-06-02.txt")

unique(gsea_complete2$comparison)
unique(gsea_complete2$timepoint)
table(gsea_complete2$timepoint, gsea_complete2$comparison)

gsea_complete2$comparison <- factor(gsea_complete2$comparison, levels=c("yo","dox"), ordered=T)

gsea_complete2$timepoint_comparison <- paste(gsea_complete2$timepoint, gsea_complete2$comparison, sep="_")

gsea_complete2$timepoint_comparison <- factor(gsea_complete2$timepoint_comparison, 
                                              levels=c("A06h_yo","A06h_dox",
                                                       "A18h_yo","A18h_dox",
                                                       "A48h_yo","A48h_dox",
                                                       "D06h_yo","D06h_dox",
                                                       "D18h_yo","D18h_dox",
                                                       "D48h_yo","D48h_dox"),
                                              ordered=T)


# Pick pathway and subset ---

PWAY <- c("GO_HISTONE_METHYLTRANSFERASE_COMPLEX","GO_CHROMATIN_ORGANIZATION", "GO_H4_HISTONE_ACETYLTRANSFERASE_COMPLEX")

data <- gsea_complete2 %>% dplyr::filter(pathway %in% PWAY) 

## Make pathway names more legible
data$pathway <- tolower(data$pathway)
data$pathway <- gsub("_", " ", data$pathway)
data$pathway <- gsub("^...", "", data$pathway)

# Dot plots
ggplot(data , aes(x=timepoint_comparison, y=pathway, size=-log(padj,10), color=NES)) + 
  geom_point(alpha = 1) + 
  scale_color_gradient2(high = "red",
                        mid = "white",
                        low = "blue") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size=10, angle=45, hjust=1), axis.title.x = element_blank(),
        axis.title.y=element_blank()) +
  theme(legend.title = element_text(size=8), 
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, 'cm'),
        legend.box = "horizontal")

ggsave(paste0("plots/yo-dox_pathways-chromatin-subset_", "timepoints_wide_dotplot.pdf"), width = 7, height = 1.5)







#========
sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.4.1
# [1] viridis_0.6.3         viridisLite_0.4.2     cowplot_1.1.1         clusterProfiler_4.8.1
# [5] NCmisc_1.2.0          msigdbr_7.5.1         org.Mm.eg.db_3.17.0   AnnotationDbi_1.62.1 
# [9] IRanges_2.34.0        S4Vectors_0.38.1      Biobase_2.60.0        BiocGenerics_0.46.0  
# [13] DT_0.28               fgsea_1.26.0          biomaRt_2.56.1        lubridate_1.9.2      
# [17] forcats_1.0.0         stringr_1.5.0         dplyr_1.1.2           purrr_1.0.1          
# [21] readr_2.1.4           tidyr_1.3.0           tibble_3.2.1          ggplot2_3.4.2        
# [25] tidyverse_2.0.0 
