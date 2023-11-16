# Model to predict chronological age of mouse based on cell type proportions in the SVZ

library(tidyverse)
library(ggthemes)
library(ggpubr)
library(cowplot)
theme_set(theme_cowplot())

setwd("~/Dropbox/10x_OSKM_2/8.Clock") #LOCAL

# Independent training data
prop_tidy <- readRDS("~/Dropbox/10x_OSKM_2/8.Clock/data/new_cell_prop.rds") #28 mice spanning 3-29 months, scRNA-seq of the SVZ (Buckley et al. 2023 Nature Aging)


#=====================
# Validate model ----
#=====================

# Leave one mouse out approach to train and evaluate model

pred_age <- c()
df <- c()

prop_tidy2 <- prop_tidy[,-12]

for(i in 1:28){
  #Remove one mouse  
  sub <- prop_tidy2[-i,] 
  
  #Make model on remaining mice
  model <- lm(data = sub, Age ~ .)
  
  #Test model on held-out mouse
  test <- predict(model, newdata = prop_tidy2[i,-11])
  
  pred_age <- c(pred_age, test)
  df <- rbind(df, prop_tidy[i,])
}

length(pred_age) #28

df$Predicted_Age <- pred_age

write.csv(df, "proportions_model_training-testing-data.csv")


# Plot
ggplot(df, aes(x = Age, y = Predicted_Age)) +
  geom_smooth(method='lm') +
  geom_point(shape=1, size=1.5) +
  ylab("Predicted age") +
  theme(axis.title.x = element_text(size=8), axis.text.x = element_text(size=8),
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8)) +
  scale_x_continuous(limits=c(0,30)) +
  stat_cor(label.y = 30, size = 2.7)  # pearson r not rho

ggsave("plots/celltype_regression_leaveonemouseout_correlation.pdf", height = 2, width = 2)



#=====================
# Model----
#=====================

model <- lm(data = prop_tidy2, Age ~ .)
summary(model)
# Call:
#   lm(formula = Age ~ ., data = prop_tidy2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -10.6892  -1.7588  -0.0121   2.5020   7.7983 
# 
# Coefficients: (1 not defined because of singularities)
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)      -94.98     214.49  -0.443    0.663
# Microglia        176.32     211.15   0.835    0.415
# Oligodendro      106.15     215.75   0.492    0.629
# Neuroblast        92.69     212.52   0.436    0.668
# Astrocyte_qNSC   132.32     226.19   0.585    0.566
# aNSC_NPC          49.74     222.85   0.223    0.826
# Endothelial       12.34     245.61   0.050    0.960
# Mural            174.44     205.06   0.851    0.406
# OPC               65.86     225.42   0.292    0.773
# Macrophage        33.24     221.89   0.150    0.883
# Ependymal            NA         NA      NA       NA
# 
# Residual standard error: 4.791 on 18 degrees of freedom
# Multiple R-squared:  0.7654,	Adjusted R-squared:  0.6481 
# F-statistic: 6.526 on 9 and 18 DF,  p-value: 0.0003814


#======================
# Apply model to whole-body OSKM data ----
#=====================

new_prop <- read.csv("data/celltypefreq_OSKM_combined.csv")
OSKM_prop <- new_prop[,c(1:11,17:19)]


# Test on data
test <- predict(model, newdata = OSKM_prop)

OSKM_prop$pred.age <- test


# Factors
OSKM_prop$Age_Treatment <- factor(OSKM_prop$Age_Treatment, 
                                  levels = c("young_untr", "old_untr", "old_2Dox0","old_2Dox5"), 
                            ordered = T)
OSKM_prop$Exp <- factor(OSKM_prop$Exp, 
                                  levels = c("1", "2"), 
                                  ordered = T)


doxagecolors <- c("turquoise2","firebrick", "mediumorchid", "plum")
names(doxagecolors) <- c("young_untr","old_untr","old_2Dox0","old_2Dox5")

write.csv(OSKM_prop, "proportions_model_OSKM-combined.csv")
OSKM_prop <- read.csv("proportions_model_OSKM-combined.csv")


# Plot
# Exclude mice not included in this study
subprop <- OSKM_prop %>% filter(Age_Treatment!="old_2Dox5")

ggplot(subprop, aes(x = Age_Treatment, y = pred.age, color=Age_Treatment)) +
  geom_point(aes(color = Age_Treatment), size = 2, shape=16, position = position_jitter(width = 0.05, height=0)) +
  scale_color_manual(values = alpha(doxagecolors,0.7), guide = 'none') +
#  labs(title = "Age prediction from celltype proportion linear regression model") +
  theme(title = element_text(size=9), axis.title.x = element_blank(), axis.title.y = element_text(size = 12), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_summary(fun = "median", geom="crossbar", width = 0.5) +
  theme(plot.title = element_text(size = 10)) +
  ylab("Predicted age") 
ggsave("plots/combined_proportions_predage-agetreatment_median.pdf", height = 2.5, width = 2.5)



#======================
# Apply model to SVZ-targeted OSKM data ----
#=====================

setwd("~/Dropbox/10x_expAC/6.Clock")
prop_tidy <- readRDS("~/Dropbox/10x_OSKM_2/8.Clock/data/new_cell_prop.rds")
new_prop <- read.csv("data/celltype_prop_SVZ-OSKM.csv")
# Mural = vascular smooth muscle cells + pericytes

# Remove 'extra' cell types (in this dataset but not in training data) and recalculate proportions ----
new_prop <- new_prop %>% mutate(Total = rowSums(new_prop[,1:10]))

OSKM_prop <- new_prop[,c(1:11,17)] %>%  mutate(Astrocyte_qNSC = Astrocyte_qNSC/Total,
                                               aNSC_NPC = aNSC_NPC/Total,
                                               Neuroblast = Neuroblast/Total,
                                               Oligodendro = Oligodendro/Total,
                                               OPC = OPC/Total,
                                               Endothelial = Endothelial/Total,
                                               Microglia = Microglia/Total,
                                               Macrophage = Macrophage/Total,
                                               Mural = Mural/Total,
                                               Ependymal = Ependymal/Total) 


# Fit and test ----
prop_tidy2 <- prop_tidy[,-12]
model <- lm(data = prop_tidy2, Age ~ .)

# Test on data ----
test <- predict(model, newdata = OSKM_prop[,1:10])

OSKM_prop$pred.age <- test

# Fix factors ----
age_treatment <- plyr::mapvalues(x = OSKM_prop$mouseID, 
                                 from = c("A.1108L.3.9.M.lane1.untr.TGTGATGG", 
                                          "B.1110L.3.9.M.lane1.untr.TCAATGGC", 
                                          "C.L0111.28.1.M.lane2.untr.CTCTAGAC", 
                                          "D.L0112.28.1.M.lane3.Dox.ACCAATGC", 
                                          "E.L0114.28.1.M.lane2.untr.AGTTGCGT", 
                                          "F.L0116.28.1.F.lane3.Dox.CGAACAAG",
                                          "G.L0117.28.1.F.lane2.untr.GTACCTGT",
                                          "H.L125.27.4.M.lane3.Dox.GAAGCTTG",
                                          "I.1114L.4.0.F.lane4.untr.AAGTACGC",
                                          "J.1115L.4.0.F.lane4.untr.ATTCGCAC",
                                          "K.1116L.4.0.F.lane4.untr.GAGTCGAT",
                                          "L.L131.27.4.F.lane6.Dox.AAGGCTAG",
                                          "M.L126.27.4.M.lane5.untr.CAGTTAGG",
                                          "N.L128.27.4.F.lane6.Dox.AACCGAAC",
                                          "O.L132.27.4.F.lane5.untr.AAGCAGTC",
                                          "P.L133.26.3.M.lane6.Dox.GAATCAGG",
                                          "Q.L0108.28.8.M.lane5.untr.ACTCGAAG"), 
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
OSKM_prop$age_treatment <- age_treatment
OSKM_prop$age_treatment <- factor(OSKM_prop$age_treatment, levels = c("young_untr", "old_untr", "old_Dox"), 
                                  ordered = T)

OSKM_prop$Exp <- c(rep("1",times=8), rep("2",times=9))
OSKM_prop$Exp <- factor(OSKM_prop$Exp, 
                        levels = c("1", "2"), 
                        ordered = T)

write.csv(OSKM_prop,"proportions_model_expAC.csv")
OSKM_prop <- read.csv("proportions_model_expAC.csv")

# Set colors ----
doxagecolors <- c("turquoise","tomato", "#612a95")
names(doxagecolors) <- unique(OSKM_prop$age_treatment)


# Plots ----

ggplot(OSKM_prop, aes(x = age_treatment, y = pred.age, color=age_treatment)) +
  geom_point(aes(color = age_treatment, shape=Exp), size = 2, position = position_jitter(width = 0.1, height=0)) +
  scale_color_manual(values = alpha(doxagecolors,0.7), guide = 'none') +
  scale_shape_manual(values=c(16,17), guide='none') +
  #  labs(title = "Age prediction from celltype proportion linear regression model") +
  theme(title = element_text(size=9), axis.title.x = element_blank(), axis.title.y = element_text(size = 12), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_summary(fun = "median", geom="crossbar", width = 0.5) +
  theme(plot.title = element_text(size = 10)) +
  ylab("Predicted age")

ggsave("plots/proportions_predage-agetreatment_shape-exp_median_small.pdf", height = 2.5, width = 2.5)



