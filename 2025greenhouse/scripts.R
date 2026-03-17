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
data <- read.csv("2025greenhousemaster_hydrogel.csv", header = TRUE, row.names = 1)
ions <- read.csv("ions/PRS_2025_hydrogel.csv", header = TRUE, row.names = 1)

# ****************************************************
# Data standardization ----
# ****************************************************
## Soil geochem ----
# Subset response and explanatory variables
treatments <- data$Treatment_Category
soil_properties <- data[ , !colnames(data) %in% "Treatment_Category"]

soil_properties_rel <- decostand(soil_properties, method = "max") # Relativize by column maxima

## Soil ions ----
# Subset response and explanatory variables
ions_treatments <- ions$Treatment_Category
ions_properties <- ions[ , !colnames(ions) %in% "Treatment_Category"]

ions_rel <- decostand(ions_properties, method = "max") # Relativize by column maxima

# ****************************************************
# Distance matrix ----
# ****************************************************
data_dist <- vegdist(soil_properties_rel, method = "euclidean")
ions_dist <- vegdist(ions_rel, method = "euclidean") 

# ****************************************************
# Graphic preparation ----
# ****************************************************
# color palette 
treatment_colors <- c(
  "BioBead" = "#CDEE40FF",  
  "Hydrogel" = "#8BD4F7FF",  
  "Microbe_Liquid" = "#E96C0CFF", 
  "Negative" = "#000000" 
)

# ****************************************************
# PCA ----
# ****************************************************
## Soil geochem ----
soil_pca <- princomp(soil_properties_rel, cor = FALSE)
summary(soil_pca, loadings = TRUE, cutoff = 0) #Proportion of Variance for axis labels below

# scree plot
screeplot(soil_pca, bstick = TRUE, type = "lines",
          main = "PCA Screeplot - Soil Properties")

pca1 <- ggbiplot(soil_pca, groups = data$Treatment_Category,
                  varname.color = "grey",   
                  arrow.color = "grey") +
  geom_point(aes(colour = data$Treatment_Category), size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = data$Treatment_Category), level = 0.95,
               linetype = "dashed", linewidth = 0.4) +
  labs(colour = "Treatment Group",
       x = "PC1 (60.4%)",   # values taken from summary()
       y = "PC2 (30.4%)",
       caption = "Figure 3: PCA of Soil Leachate and Moisture") +
  theme_bw()

pca1

## Soil ions ----
ion_pca <- prcomp(ions_rel) #can't do prcomp with fewer units than variables
summary(ion_pca, loadings = TRUE, cutoff = 0) #Proportion of Variance for axis labels below
ion_pca$rotation #for reports on loading

# scree plot
screeplot(ion_pca, bstick = TRUE, type = "lines",
          main = "PCA Screeplot - Soil Ions")

pca2 <- ggbiplot(ion_pca, groups = ions$Treatment_Category,
                 varname.color = "grey",   
                 arrow.color = "grey") +
  geom_point(aes(colour = ions$Treatment_Category), size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = ions$Treatment_Category), level = 0.95,
               linetype = "dashed", linewidth = 0.4) +
  labs(colour = "Treatment Group",
       x = "PC1 (36.4%)",   # values taken from summary()
       y = "PC2 (21.1%)",
       caption = "Figure 4: PCA of Soil Ion Supply Rate") +
  theme_bw()

pca2

# ****************************************************
# NMDS ----
# ****************************************************
## Soil geochem ----
set.seed(1)
nmds_result <- metaMDS(soil_properties_rel, 
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
nmds1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = cluster), level = 0.95, linetype = "dashed", linewidth = 0.4) +  # potentially correct to ggforce (below), this assumes normal dist
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "NMDS1", y = "NMDS2",
       caption = "Figure 5: NMDS Ordination of Soil Leachate and Moisture")

nmds1

## Soil ions ----
set.seed(1)
nmds_ions <- metaMDS(ions_dist,
                     k = 2,
                     trymax = 100)
nmds_ions$stress

#create data frame
cluster_df_ions <- data.frame( 
  label = rownames(ions),
  cluster = as.factor(ions$Treatment_Category)
)
head(cluster_df_ions)

nmds_scores_ions <- vegan::scores(nmds_ions, tidy = TRUE) |> 
  filter(score == "sites") |> 
  left_join(cluster_df_ions, by = "label")

nmds2 <- ggplot(nmds_scores_ions, aes(x = NMDS1, y = NMDS2, colour = cluster)) +
  geom_point(size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = cluster), level = 0.95, linetype = "dashed", linewidth = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(x = "NMDS1", y = "NMDS2",
       caption = "Figure 6: NMDS Ordination of Soil Ion Supply Rates")
nmds2

# ****************************************************
# One-way PERMANOVA ----
# ****************************************************
## Soil geochem ----
set.seed(1)
gel_permanova_result <- adonis2(data_dist ~ treatments, 
                                permutations = 9999)
gel_permanova_result

## Soil ions ----
set.seed(1)
ion_permanova_result <- adonis2(ions_dist ~ ions_treatments, 
                                permutations = 9999)
ion_permanova_result

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

## Soil ions ----
ion_gel2factor <- read.csv("ions/PRS_2025_gel2factor.csv", header = TRUE, row.names = 1)

# Subset 
ion_properties <- ion_gel2factor[ , !colnames(ion_gel2factor) %in% 
                                    c("Treatment_Category", "Hydrogel", "Microbe")]

# Relativize/redo distance matrix
ion_properties_rel <- decostand(ion_properties, method = "max")
ion2factor_dist <- vegdist(ion_properties_rel, method = "euclidean")

set.seed(1)
ion_permanova_2way <- adonis2(
  ion2factor_dist ~ Hydrogel * Microbe,
  data = ion_gel2factor,
  permutations = 9999,
  by = "term"
)
ion_permanova_2way

# ****************************************************
# RDA ----
# ****************************************************
## Soil geochem ----
gel_rda <- rda(
  soil_properties_rel ~ Hydrogel * Microbe,
  data = data_gel2factor,
  scale = FALSE  
)

summary(gel_rda)

# Significance test of the full model
anova(gel_rda, permutations = 9999)

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
anova(ion_rda, permutations = 9999)
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