# Manhattan-like plots highlighting DEGs by timepoint using DESeq2 DE values.

rm(list = ls())
library(tidyverse)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(scales)
library(NCmisc) 
theme_set(theme_cowplot())

# See sessionInfo() at end of script for version information.

# Modify path below. All subsequent paths are relative.
setwd("~/Dropbox/bulkRNAseq_expAL/1.DESeq2")


#==================================================================================================
# old untr vs. young untr ----
#==================================================================================================

df1 <- read_csv("data/DESeqRes_designCondition_A18h-old-untr_A18h-young-untr.csv")
df1$Timepoint <- "A18h"

df2 <- read_csv("data/DESeqRes_designCondition_D06h-old-untr_D06h-young-untr.csv")
df2$Timepoint <- "D06h"

df3 <- read_csv("data/DESeqRes_designCondition_D18h-old-untr_D18h-young-untr.csv")
df3$Timepoint <- "D18h"

df4 <- read_csv("data/DESeqRes_designCondition_D48h-old-untr_D48h-young-untr.csv")
df4$Timepoint <- "D48h"

df_old_young <- rbind(df1,df2,df3,df4)


# Remove NAs
df_old_young <- df_old_young[is.na(df_old_young$padj)==FALSE,]

# Add z-score based on two sided null hypothesis.
df_old_young$z <- p.to.Z(df_old_young$pvalue) * sign(df_old_young$log2FoldChange)

# Randomize rows to reduce overplotting issues.
df_old_young <- df_old_young[sample(nrow(df_old_young)), ]

# Adjust factors
df_old_young$Timepoint <- factor(df_old_young$Timepoint)

# Adjust colors
#timecols <- inferno(6)
timecols <- c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF" )
names(timecols) <- c("A06h", "A18h", "A48h", "D06h", "D18h", "D48h")


#================
# P_ADJ < 0.1

# Make custom color column to facilitate grey coloring by threshold.
col <- timecols[as.character(df_old_young$Timepoint)]
col[df_old_young$padj > 0.1] <- "#D3D3D3" # grey
df_old_young$col <- as.factor(col)

q <- ggplot(df_old_young, aes(x = Timepoint, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .4, size = 1, shape=16) +
  #scale_shape_manual(values=rep(16,nlevels(df_old_young$Timepoint))) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Old vs. Young") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df_old_young$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q

png('plots_manhattan/manhattan_old.young_padj-0.1.png', width=5, height=5, units = 'in', res=300)
q
dev.off()


#==================================================================================================
# old Dox vs. old untr ----
#==================================================================================================

df1 <- read_csv("data/DESeqRes_designCondition_A18h-old-Dox_A18h-old-untr.csv")
df1$Timepoint <- "A18h"

df2 <- read_csv("data/DESeqRes_designCondition_D06h-old-Dox_D06h-old-untr.csv")
df2$Timepoint <- "D06h"

df3 <- read_csv("data/DESeqRes_designCondition_D18h-old-Dox_D18h-old-untr.csv")
df3$Timepoint <- "D18h"

df4 <- read_csv("data/DESeqRes_designCondition_D48h-old-Dox_D48h-old-untr.csv")
df4$Timepoint <- "D48h"

df_dox_untr <- rbind(df1,df2,df3,df4)


# Remove NAs
df_dox_untr <- df_dox_untr[is.na(df_dox_untr$padj)==FALSE,]

# Add z-score based on two sided null hypothesis.
df_dox_untr$z <- p.to.Z(df_dox_untr$pvalue) * sign(df_dox_untr$log2FoldChange)
df_dox_untr$z.adj <- p.to.Z(df_dox_untr$padj) * sign(df_dox_untr$log2FoldChange)

# Randomize rows to reduce overplotting issues.
df_dox_untr <- df_dox_untr[sample(nrow(df_dox_untr)), ]

# Adjust factors
df_dox_untr$Timepoint <- factor(df_dox_untr$Timepoint)

# Adjust colors
#timecols <- inferno(6)
timecols <- c("#000004FF", "#330A5FFF", "#781C6DFF", "#BB3754FF", "#ED6925FF", "#FCB519FF" )
names(timecols) <- c("A06h", "A18h", "A48h", "D06h", "D18h", "D48h")

#================
# P_ADJ < 0.1

# Adjust colors

# Make custom color column to facilitate grey coloring by threshold.
col <- timecols[as.character(df_dox_untr$Timepoint)]
col[df_dox_untr$padj > 0.1] <- "#D3D3D3" # grey
df_dox_untr$col <- as.factor(col)

q <- ggplot(df_dox_untr, aes(x = Timepoint, y = z, color = col)) +
  geom_jitter(width = 0.40, alpha = .5, size = 1, shape=16) +
  #scale_shape_manual(values=rep(16,nlevels(df_dox_untr$Timepoint))) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 14), axis.title.x = element_blank()) +
  ggtitle("Dox vs. untr") +
  theme(axis.title.y = element_text(size = 20, face = "plain")) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 18)) +
  theme(plot.title = element_text(size=20, face = "plain")) +
  labs(y = "Z-score") +
  theme(legend.position="none") +
  scale_color_manual(values = levels(df_dox_untr$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q

png('plots_manhattan/manhattan_dox.untr_padj-0.1.png', width=5, height=5, units = 'in', res=300)
q
dev.off()



#=============
sessionInfo()
# [1] NCmisc_1.2.0       scales_1.2.1       viridis_0.6.3      viridisLite_0.4.2  RColorBrewer_1.1-3
# [6] cowplot_1.1.1      lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.2       
# [11] purrr_1.0.1        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.2     
# [16] tidyverse_2.0.0  