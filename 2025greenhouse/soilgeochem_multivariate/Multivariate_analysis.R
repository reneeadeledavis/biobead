# ****************************************************
# Load packages ----
# ****************************************************
library(vegan)
library(dplyr)
library(tidyr)
library(ggforce) #soil ion plots/ellipses
library(tidyverse)
library(ggbiplot)
library(ggplot2)

# ****************************************************
# Read in data ----
# ****************************************************
data <- read.csv("2025greenhousemaster_alleffects.csv", header = TRUE, row.names = 1)

# ****************************************************
# Data standardization ----
# ****************************************************
# Subset response and explanatory variables
treatments <- data$Treatment_Category
properties <- data[ , !colnames(data) %in% "Treatment_Category"]

properties_rel <- decostand(properties, method = "max") # Relativize by column maxima

# ****************************************************
# Distance matrix ----
# ****************************************************
data_dist <- vegdist(properties_rel, method = "euclidean")

# ****************************************************
# Graphic preparation ----
# ****************************************************
# color palette 
treatment_colors <- c(
  "BioBead"        = "#8DD3C7",
  "BioBead_SL"     = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Liquid_Microbes"= "#FFFFB3",  
  "HardyGro"       = "#FDB462",
  "Hydrogel"       = "#80B1D3",
  "Negative"       = "#D9D9D9"    
)

# ****************************************************
# PCA ----
# ****************************************************
pca <- princomp(properties_rel, cor = FALSE)
summary(pca, loadings = TRUE, cutoff = 0) #Proportion of Variance for axis labels below

# scree plot
screeplot(pca, bstick = TRUE, type = "lines",
          main = "PCA Screeplot") #retain 4 components

pca_plot <- ggbiplot(pca, groups = data$Treatment_Category,
                 varname.color = "grey",   
                 arrow.color = "grey") +
  geom_point(aes(colour = data$Treatment_Category), size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = data$Treatment_Category), level = 0.95,
               linetype = "dashed", linewidth = 0.4) +
  labs(colour = "Treatment Group",
       x = "PC1 (33.3%)",   # values taken from summary()
       y = "PC2 (25.0%)",
       caption = "Figure 1: PCA") +
  theme_bw()

pca_plot

# ****************************************************
# NMDS ----
# ****************************************************
set.seed(1)
nmds_result <- metaMDS(properties_rel, 
                       distance = "euclidean", 
                       k = 3, 
                       trymax = 100)
nmds_result$stress

#create data frame
cluster_df <- data.frame( 
  label = rownames(data),
  cluster = as.factor(data$Treatment_Category)
)
head(cluster_df)

nmds_scores <- vegan::scores(nmds_result, tidy = TRUE) |> 
  filter(score == "sites") |> 
  left_join(cluster_df, by = "label")

#plot
nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = cluster), level = 0.95, linetype = "dashed", linewidth = 0.4) +  # potentially correct to ggforce (below), this assumes normal dist
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "NMDS1", y = "NMDS2",
       caption = "Figure 2: NMDS Ordination")

nmds

# ****************************************************
# One-way PERMANOVA ----
# ****************************************************
## FIX ALL THIS ----
set.seed(1)
gel_permanova_result <- adonis2(data_dist ~ treatments, 
                                permutations = 9999)
gel_permanova_result


# ****************************************************
# Two-factor PERMANOVA ----
# ****************************************************
## Soil geochem ----
data_gel2factor <- read.csv("2025greenhousemaster_gel2factor.csv", header = TRUE, row.names = 1)

# Subset 
soil_properties <- data_gel2factor[ , !colnames(data_gel2factor) %in% c("Treatment_Category", "Hydrogel", "Microbe")]

# Relativize/redo distance matrix
soil_properties_rel <- decostand(soil_properties, method = "max")
gel2factor_dist <- vegdist(soil_properties_rel, method = "euclidean")

set.seed(1)
gel_permanova_2way <- adonis2(
  gel2factor_dist ~ Hydrogel * Microbe,
  data = data_gel2factor,
  permutations = 9999,
  by = "term"
)
gel_permanova_2way


# ****************************************************
# RDA ----
# ****************************************************
## Soil geochem ----
gel_rda <- rda(
  properties_rel ~ Hydrogel * Microbe,
  data = data_gel2factor,
  scale = FALSE  
)

summary(gel_rda)

# Significance test of the full model
set.seed(1)
anova(gel_rda, by = "term", permutations = 9999) 

# Plot
plot(gel_rda, type = "n")
points(gel_rda, display = "sites", 
       col = treatment_colors[data_gel2factor$Treatment_Category], 
       pch = 16, cex = 1.2)
text(gel_rda, display = "bp", col = "#1927DDFF")  # arrows for explanatory variables
legend("topright", legend = names(treatment_colors), 
       col = treatment_colors, pch = 16, cex = 0.7)
mtext("Figure 7: RDA of Soil Leachate and Moisture", 
      side = 1, # bottom 
      line = 4, # how far below the axis
      cex = 0.8) # font size

## Soil ions ----
ion_rda <- rda(
  ion_properties_rel ~ Hydrogel * Microbe,
  data = ion_gel2factor,
  scale = FALSE  
)
summary(ion_rda)

# Significance test of the full model
set.seed(1)
anova(ion_rda, by = "term", permutations = 9999)

# Plot (run together)
plot(ion_rda, type = "n")
points(ion_rda, display = "sites", 
       col = treatment_colors[ion_gel2factor$Treatment_Category], 
       pch = 16, cex = 1.2)
text(ion_rda, display = "bp", col = "#1927DDFF")  # arrows for explanatory variables
legend("bottomleft", legend = names(treatment_colors), 
       col = treatment_colors, pch = 16, cex = 0.7)
mtext("Figure 8: RDA of Soil Ion Supply Rate constrained by Hydrogel and Microbe presence", 
      side = 1, # bottom 
      line = 4, # how far below the axis
      cex = 0.8) # font size