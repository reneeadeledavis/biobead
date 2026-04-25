# ****************************************************
# Load packages ----
# ****************************************************
library(ggplot2)
library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)
library(vegan)
library(ggbiplot)

# ****************************************************
# Load data ----
# ****************************************************
data <- read.csv("publication/2025greenhousemaster.csv")

# ****************************************************
# Summary statistics ----
# ****************************************************
summary_data <- data %>%
  group_by(Treatment) %>%
  summarise(
    
    # GrainWeight
    mean_grain_weight = mean(GrainWeight_g, na.rm = TRUE),
    sd_grain_weight = sd(GrainWeight_g, na.rm = TRUE),
    se_grain_weight = sd(GrainWeight_g, na.rm = TRUE) / sqrt(n()),
    
    # Harvest_Index
    mean_harvest_index = mean(Harvest_Index, na.rm = TRUE),
    sd_harvest_index = sd(Harvest_Index, na.rm = TRUE),
    se_harvest_index = sd(Harvest_Index, na.rm = TRUE) / sqrt(n()),
    
    # Root_biomass
    mean_root_biomass = mean(Root_biomass_g, na.rm = TRUE),
    sd_root_biomass = sd(Root_biomass_g, na.rm = TRUE),
    se_root_biomass = sd(Root_biomass_g, na.rm = TRUE) / sqrt(n()),
    
    # NH4-N
    mean_nh4 = mean(NH4.N_ugml, na.rm = TRUE),
    sd_nh4 = sd(NH4.N_ugml, na.rm = TRUE),
    se_nh4 = sd(NH4.N_ugml, na.rm = TRUE) / sqrt(n()),
    
    # NO3-N
    mean_no3 = mean(NO3.N_ugml, na.rm = TRUE),
    sd_no3 = sd(NO3.N_ugml, na.rm = TRUE),
    se_no3 = sd(NO3.N_ugml, na.rm = TRUE) / sqrt(n()),
    
    # PO4-P
    mean_po4 = mean(PO4.P_ugml, na.rm = TRUE),
    sd_po4 = sd(PO4.P_ugml, na.rm = TRUE),
    se_po4 = sd(PO4.P_ugml, na.rm = TRUE) / sqrt(n()),
    
    # AMF colonization
    mean_pct_total = mean(AMF_colonization, na.rm = TRUE),
    sd_pct_total   = sd(AMF_colonization, na.rm = TRUE),
    se_pct_total   = sd(AMF_colonization, na.rm = TRUE) / sqrt(n()),
    
    mean_pct_arb   = mean(Arb_colonization, na.rm = TRUE),
    sd_pct_arb     = sd(Arb_colonization, na.rm = TRUE),
    se_pct_arb     = sd(Arb_colonization, na.rm = TRUE) / sqrt(n()),
    
    mean_pct_hyph  = mean(Hyph_colonization, na.rm = TRUE),
    sd_pct_hyph    = sd(Hyph_colonization, na.rm = TRUE),
    se_pct_hyph    = sd(Hyph_colonization, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

# Display summary statistics
cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)

# Set color palette
treatment_colors <- c(
  "BB_AMF"         = "#8DD3C7",
  "BB_Bac"         = "#BC80BD",
  "BioBead"        = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Microbes_Liquid"= "#FFFFB3",  
  "HardyGro"       = "#FDB462",
  "Hydrogel"       = "#80B1D3",
  "Neg_ctrl"       = "#D9D9D9"    
)


# ****************************************************
# Part 1: Wheat yield ----
# ****************************************************

## Grain weight ----
### ANOVA ----
anova_grainweight <- aov(GrainWeight_g ~ Treatment, data = data)
summary(anova_grainweight)

tukey_grainweight <- TukeyHSD(anova_grainweight)

# Shapiro Wilk test for normality
shapiro.test(anova_grainweight$residuals)
# W = 0.96388, p-value = 0.09176 - data normally distributed

# Residuals distribution
hist(anova_grainweight$residuals, 
     main = "Grain Weight Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_grainweight$residuals, main = "Normal Q-Q Plot - Grain Weight Residuals")
qqline(anova_grainweight$residuals, col = "lightblue", lwd = 2)

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_grain_weight), 
                         y = mean_grain_weight,
                         fill = Treatment)) +  
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +  
  scale_fill_manual(values = treatment_colors) +  
  geom_errorbar(aes(ymin = mean_grain_weight - se_grain_weight, 
                    ymax = mean_grain_weight + se_grain_weight),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("Grain Weight")),
    y = expression(bold("Grain Weight (g)")),
    x = NULL
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

## Harvest Index ----
### ANOVA ----
anova_harvestindex <- aov(Harvest_Index ~ Treatment, data = data)
summary(anova_harvestindex)

tukey_harvestindex <- TukeyHSD(anova_harvestindex)
tukey_harvestindex

# Shapiro Wilk test for normality
shapiro.test(anova_harvestindex$residuals)
# W = 0.97042, p-value = 0.1837 - normal data

# Residuals distribution
hist(anova_harvestindex$residuals, 
     main = "Harvest Index Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_harvestindex$residuals, main = "Normal Q-Q Plot - Harvest Index Residuals")
qqline(anova_harvestindex$residuals, col = "lightblue", lwd = 2)

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_harvest_index), 
                         y = mean_harvest_index,
                         fill = Treatment)) +  
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +  
  scale_fill_manual(values = treatment_colors) +  
  geom_errorbar(aes(ymin = mean_harvest_index - se_harvest_index, 
                    ymax = mean_harvest_index + se_harvest_index),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("Harvest Index")),
    y = expression(bold("Harvest Index")),
    x = NULL
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

## Leachate ----
### ANOVA ----
#### NH4 -----
anova_nh4 <- aov(NH4.N_ugml ~ Treatment, data = data)
summary(anova_nh4)

tukey_nh4 <- TukeyHSD(anova_nh4)
tukey_nh4

# Shapiro-Wilk test for normality
shapiro.test(anova_nh4$residuals)
# W = 0.73369, p-value = 9.594e-09 - data not normal

# Residuals distribution
hist(anova_nh4$residuals,
     main = "NH4-N Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plot
qqnorm(anova_nh4$residuals, main = "Normal Q-Q Plot - NH4-N Residuals")
qqline(anova_nh4$residuals, col = "lightblue", lwd = 2)

#### NO3 ----
anova_no3 <- aov(NO3.N_ugml ~ Treatment, data = data)
summary(anova_no3)

tukey_no3 <- TukeyHSD(anova_no3)
tukey_no3

# Shapiro-Wilk test for normality
shapiro.test(anova_no3$residuals)
# W = 0.74005, p-value = 1.294e-08 - data not normal

# Residuals distribution
hist(anova_no3$residuals,
     main = "NO3-N Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plot
qqnorm(anova_no3$residuals, main = "Normal Q-Q Plot - NO3-N Residuals")
qqline(anova_no3$residuals, col = "lightblue", lwd = 2)

#### PO4 -----
anova_po4 <- aov(PO4.P_ugml ~ Treatment, data = data)
summary(anova_po4)

tukey_po4 <- TukeyHSD(anova_po4)
tukey_po4

# Shapiro-Wilk test for normality
shapiro.test(anova_po4$residuals)
# W = 0.85173, p-value = 6.446e-06 - data not normal

# Residuals distribution
hist(anova_po4$residuals,
     main = "PO4-P Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plot
qqnorm(anova_po4$residuals, main = "Normal Q-Q Plot - PO4-P Residuals")
qqline(anova_po4$residuals, col = "lightblue", lwd = 2)

### Faceted bar plot ----
leachate_long <- summary_data %>%
  pivot_longer(
    cols      = c(mean_nh4, mean_no3, mean_po4),
    names_to  = "Nutrient",
    values_to = "Mean_ugml"
  ) %>%
  mutate(
    SE = case_when(
      Nutrient == "mean_nh4" ~ se_nh4,
      Nutrient == "mean_no3" ~ se_no3,
      Nutrient == "mean_po4" ~ se_po4
    ),
    Nutrient = recode(Nutrient,
                      "mean_nh4" = "NH4-N",
                      "mean_no3" = "NO3-N",
                      "mean_po4" = "PO4-P"
    )
  )

ggplot(leachate_long, aes(x = Treatment, y = Mean_ugml, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean_ugml - SE, ymax = Mean_ugml + SE),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  facet_wrap(~ Nutrient, scales = "free_y") +
  scale_fill_manual(values = treatment_colors) +
  labs(
    title = expression(bold("NH"[4]*"-N, NO"[3]*"-N, and PO"[4]*"-P Leachate by Treatment")),
    x     = NULL,
    y     = expression(bold("Mean Concentration (µg/mL)")),
    fill  = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y       = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title         = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "none"
  )

# ****************************************************
# Part 2: AMF establishment ----
# ****************************************************

## Stacked bar plot ----
summary_stacked <- data %>%
  group_by(Treatment) %>%
  summarise(
    mean_Arb   = mean(Arb_colonization,  na.rm = TRUE),
    mean_Hyph  = mean(Hyph_colonization, na.rm = TRUE),
    se_total   = sd(AMF_colonization, na.rm = TRUE) / sqrt(sum(!is.na(AMF_colonization))),
    total_mean = mean(AMF_colonization,  na.rm = TRUE),
    .groups = "drop"
  )

stacked_long <- summary_stacked %>%
  pivot_longer(
    cols      = c(mean_Arb, mean_Hyph),
    names_to  = "Structure",
    values_to = "Mean_Pct"
  ) %>%
  mutate(
    Structure = recode(Structure,
                       "mean_Arb"  = "Arbuscule",
                       "mean_Hyph" = "Intraradical Hyphae"
    ),
    Structure = factor(Structure,
                       levels = c("Intraradical Hyphae", "Arbuscule")
    )
  )

structure_colors <- c(
  "Intraradical Hyphae" = "#80B1D3",
  "Arbuscule"           = "#BEBADA"
)

ggplot(stacked_long, aes(x = Treatment, y = Mean_Pct, fill = Structure)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(
    data = summary_stacked,
    aes(x = Treatment, ymin = total_mean - se_total, ymax = total_mean + se_total),
    inherit.aes = FALSE,
    width = 0.2, linewidth = 0.6, color = "grey30"
  ) +
  scale_fill_manual(values = structure_colors) +
  labs(
    title = expression(bold("AMF Root Colonization by Structure")),
    x     = NULL,
    y     = expression(bold("Mean Percent Colonization (%)")),
    fill  = "AMF Structure"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y       = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title         = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

## Total colonization ----
### ANOVA ----
anova_colonization <- aov(AMF_colonization ~ Treatment,
                          data = data)
summary(anova_colonization)

tukey_colonization <- TukeyHSD(anova_colonization)
print(tukey_colonization)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_colonization$residuals)
# W = 0.94717, p-value = 0.01593 - data not normal

# Residuals histogram
hist(anova_colonization$residuals,
     main = "Percent Total Colonization Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_colonization$residuals, main = "Normal Q-Q Plot - Colonization Residuals")
qqline(anova_colonization$residuals, col = "lightblue", lwd = 2)

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_pct_total),
                         y = mean_pct_total,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_total - se_pct_total,
                    ymax = mean_pct_total + se_pct_total),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("AMF Total Root Colonization")),
    x = NULL,
    y = expression(bold("Percent Colonization (%)"))
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

## Arbuscule establishment ----
### ANOVA ----
anova_arbuscule <- aov(Arb_colonization ~ Treatment,
                       data = data)
summary(anova_arbuscule)

tukey_arbuscule <- TukeyHSD(anova_arbuscule)
print(tukey_arbuscule)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_arbuscule$residuals)
# W = 0.9693, p-value = 0.1632, data normal

# Residuals histogram
hist(anova_arbuscule$residuals,
     main = "Arbuscule Colonization Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_arbuscule$residuals, main = "Normal Q-Q Plot - Colonization Residuals")
qqline(anova_arbuscule$residuals, col = "lightblue", lwd = 2)

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_pct_arb),
                         y = mean_pct_arb,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_arb - se_pct_arb,
                    ymax = mean_pct_arb + se_pct_arb),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("Arbuscule Establishment")),
    x = NULL,
    y = expression(bold("Arbuscule Colonization (%)"))
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

## Intra hyphae establishment ----
### ANOVA ----
anova_hyph <- aov(Hyph_colonization ~ Treatment,
                  data = data)
summary(anova_hyph)

tukey_hyph <- TukeyHSD(anova_hyph)
print(tukey_hyph)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_hyph$residuals)
# W = 0.95547, p-value = 0.03757- data borderline normal

# Residuals histogram
hist(anova_hyph$residuals,
     main = "Intraradical Hyphae Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_hyph$residuals, main = "Normal Q-Q Plot - Intraradical Hyphae Residuals")
qqline(anova_hyph$residuals, col = "lightblue", lwd = 2)

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_pct_hyph),
                         y = mean_pct_hyph,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_hyph - se_pct_hyph,
                    ymax = mean_pct_hyph + se_pct_hyph),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("Intraradical Hyphae Establishment")),
    x = NULL,
    y = expression(bold("Intraradical Hyphae (%)"))
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ****************************************************
# Part 3: Correlations ----
# ****************************************************

## Pearson's correlation ----
cor_vars <- data[, c("AMF_colonization", "Arb_colonization", "Hyph_colonization",
                       "GrainWeight_g", "Root_biomass_g", "Total_biomass")]

cor(cor_vars, use = "complete.obs", method = "pearson")

## test
result <- rcorr(as.matrix(cor_vars), type = "pearson")
result$r   # correlation coefficients
result$P   # p-values

## plot
corrplot(result$r,
         method   = "color",
         type     = "upper",
         tl.col   = "black",
         tl.srt   = 45,
         addCoef.col = "black",
         number.cex  = 0.8,
         p.mat    = result$P,
         sig.level = 0.05,
         insig    = "blank")

## PCA ----
## data standardization
treatments <- data$Treatment

properties <- data[, c("AMF_colonization", "Arb_colonization", "Hyph_colonization",
                         "GrainWeight_g", "Root_biomass_g", "Total_biomass")]

properties_rel <- decostand(properties, method = "max")

## PCA
pca <- princomp(properties_rel, cor = FALSE)
summary(pca, loadings = TRUE, cutoff = 0)

# Scree plot
screeplot(pca, bstick = TRUE, type = "lines",
          main = "PCA Screeplot")

# Biplot — update axis label percentages from summary() output
pca_plot <- ggbiplot(pca, groups = treatments,
                     varname.color = "grey40",
                     arrow.color   = "grey40") +
  geom_point(aes(colour = treatments), size = 3) +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  stat_ellipse(aes(group = treatments), level = 0.95,
               linetype = "dashed", linewidth = 0.4) +
  labs(colour = "Treatment",
       x = "PC1 (66.7%)",   # update from summary()
       y = "PC2 (13.4%)",
       caption = "PCA: AMF colonization and plant response variables") +
  theme_bw()

pca_plot

# focus on biological groups
group_map <- c(
  "BB_AMF"         = "AMF Beads",
  "BioBead"       = "AMF Beads",
  "BB_Bac"         = "Non-AMF Beads",
  "Hydrogel"          = "Non-AMF Beads",
  "Microbes_Liquid"        = "Non-Bead Controls",
  "HardyGro" = "Non-Bead Controls",
  "Fertilizer"          = "Non-Bead Controls",
  "Neg_ctrl"            = "Non-Bead Controls"
)

groups_collapsed <- group_map[treatments]

group_colors <- c(
  "AMF Beads"         = "#BEBADA",
  "Non-AMF Beads"     = "#8DD3C7",
  "Non-Bead Controls" = "#D9D9D9"
)

pca_plot_groups <- ggbiplot(pca, groups = groups_collapsed,
                            varname.color = "grey40",
                            arrow.color   = "grey40") +
  geom_point(aes(colour = groups_collapsed), size = 3) +
  scale_color_manual(values = group_colors, name = "Treatment Group") +
  stat_ellipse(aes(group = groups_collapsed), level = 0.95,
               linetype = "dashed", linewidth = 0.4) +
  labs(colour = "Treatment Group",
       x = "PC1 (66.7%)",
       y = "PC2 (13.4%)",
       caption = "PCA: AMF colonization and plant response variables") +
  theme_bw()

pca_plot_groups
