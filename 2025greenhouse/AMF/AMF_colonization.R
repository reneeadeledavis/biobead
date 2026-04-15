# ****************************************************
# Load packages ----
# ****************************************************
library(ggstatsplot)
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(dplyr)
library(plotrix)
library(see)
library(ggpubr)

# ****************************************************
# Read in data ----
# ****************************************************

colonization_2025 <- read.csv("AMF/Root_colonization_2025.csv", header = TRUE) %>%
  filter(!is.na(Treatment_category) & Treatment_category != "")

# ****************************************************
# Summary statistics ----
# ****************************************************

## complete data set
summary_data_complete <- colonization_2025 %>%
  group_by(Treatment_category) %>%
  summarise(
    
    # Total colonization
    mean_pct_total = mean(Percent_colonization_total, na.rm = TRUE),
    sd_pct_total   = sd(Percent_colonization_total, na.rm = TRUE),
    se_pct_total   = sd(Percent_colonization_total, na.rm = TRUE) / sqrt(n()),
    
    # Vesicle/Spore colonization
    mean_pct_vs    = mean(Percent_colonization_VS, na.rm = TRUE),
    sd_pct_vs      = sd(Percent_colonization_VS, na.rm = TRUE),
    se_pct_vs      = sd(Percent_colonization_VS, na.rm = TRUE) / sqrt(n()),
    
    # Arbuscule colonization
    mean_pct_arb   = mean(Percent_colonization_Arb, na.rm = TRUE),
    sd_pct_arb     = sd(Percent_colonization_Arb, na.rm = TRUE),
    se_pct_arb     = sd(Percent_colonization_Arb, na.rm = TRUE) / sqrt(n()),
    
    # Intraradical hyphae colonization
    mean_pct_hyph  = mean(Percent_colonization_hyphae, na.rm = TRUE),
    sd_pct_hyph    = sd(Percent_colonization_hyphae, na.rm = TRUE),
    se_pct_hyph    = sd(Percent_colonization_hyphae, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

cat("=== Summary Statistics ===\n")
print(summary_data_complete, n = Inf)

## filtered data set
summary_data <- colonization_filtered %>%
  group_by(Treatment_category) %>%
  summarise(
    
    # Total colonization
    mean_pct_total = mean(Percent_colonization_total, na.rm = TRUE),
    sd_pct_total   = sd(Percent_colonization_total, na.rm = TRUE),
    se_pct_total   = sd(Percent_colonization_total, na.rm = TRUE) / sqrt(n()),
    
    # Vesicle/Spore colonization
    mean_pct_vs    = mean(Percent_colonization_VS, na.rm = TRUE),
    sd_pct_vs      = sd(Percent_colonization_VS, na.rm = TRUE),
    se_pct_vs      = sd(Percent_colonization_VS, na.rm = TRUE) / sqrt(n()),
    
    # Arbuscule colonization
    mean_pct_arb   = mean(Percent_colonization_Arb, na.rm = TRUE),
    sd_pct_arb     = sd(Percent_colonization_Arb, na.rm = TRUE),
    se_pct_arb     = sd(Percent_colonization_Arb, na.rm = TRUE) / sqrt(n()),
    
    # Intraradical hyphae colonization
    mean_pct_hyph  = mean(Percent_colonization_hyphae, na.rm = TRUE),
    sd_pct_hyph    = sd(Percent_colonization_hyphae, na.rm = TRUE),
    se_pct_hyph    = sd(Percent_colonization_hyphae, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)

# ****************************************************
# Color palette ----
# ****************************************************

##complete data set
treatment_colors_complete <- c(
  "BB_Bacteria"        = "#8DD3C7",
  "BB_AMF"             = "#BEBADA",
  "BB_Mix"             = "#FB8072",
  "SL_BB_Mix"          = "#FFFFB3",
  "SL_BB_Bac"          = "#FDB462",
  "SL_BB_AMF"          = "#80B1D3",
  "BB_Mix_Reduced"     = "#B3DE69",
  "HalfFert_HalfBB"    = "#FCCDE5",
  "Liquid_consortium"  = "#D9D9D9",
  "HardyGro"           = "#BC80BD",
  "Empty_Bead"         = "#CCEBC5",
  "SL_Empty_Bead"      = "#FFED6F",
  "Fertilizer_control" = "#FC8D62",
  "Negative_control"   = "#E5C494"
)

##focused data set
treatment_colors <- c(
  "BB_AMF"    = "#8DD3C7",
  "BB_Bac"    = "#BC80BD",
  "BioBead"        = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Microbes_Liquid"= "#FFFFB3",  
  "HardyGro"       = "#FDB462",
  "Hydrogel"       = "#80B1D3",
  "Neg_ctrl"       = "#D9D9D9"    
)

# ****************************************************
# Percent colonization bar plot ----
# ****************************************************

## filtered data set
ggplot(summary_data, aes(x = reorder(Treatment_category, mean_pct_total),
                         y = mean_pct_total,
                         fill = Treatment_category)) +
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

# ****************************************************
# Percent arbuscule bar plot ----
# ****************************************************

## filtered data set
ggplot(summary_data, aes(x = reorder(Treatment_category, mean_pct_arb),
                         y = mean_pct_arb,
                         fill = Treatment_category)) +
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

# ****************************************************
# Percent intra hyphae bar plot ----
# ****************************************************
ggplot(summary_data, aes(x = reorder(Treatment_category, mean_pct_hyph),
                         y = mean_pct_hyph,
                         fill = Treatment_category)) +
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
# ANOVA ----
# ****************************************************

## Arbuscule % colonization
anova_arbuscule <- aov(Percent_colonization_Arb ~ Treatment_category,
                          data = colonization_filtered)
summary(anova_arbuscule)

tukey_arbuscule <- TukeyHSD(anova_arbuscule)
print(tukey_arbuscule)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_arbuscule$residuals)

 # Residuals histogram
hist(anova_arbuscule$residuals,
     main = "Arbuscule Colonization Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_arbuscule$residuals, main = "Normal Q-Q Plot - Colonization Residuals")
qqline(anova_arbuscule$residuals, col = "lightblue", lwd = 2)

## Total % colonization
anova_colonization <- aov(Percent_colonization_total ~ Treatment_category,
                          data = colonization_filtered)
summary(anova_colonization)

tukey_colonization <- TukeyHSD(anova_colonization)
print(tukey_colonization)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_colonization$residuals)

# Residuals histogram
hist(anova_colonization$residuals,
     main = "Percent Colonization Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_colonization$residuals, main = "Normal Q-Q Plot - Colonization Residuals")
qqline(anova_colonization$residuals, col = "lightblue", lwd = 2)

## Intraradical hyphae establishment
anova_hyph <- aov(Percent_colonization_hyphae ~ Treatment_category,
                       data = colonization_filtered)
summary(anova_hyph)

tukey_hyph <- TukeyHSD(anova_hyph)
print(tukey_hyph)

# Shapiro-Wilk test for normality of residuals
shapiro.test(anova_hyph$residuals)

# Residuals histogram
hist(anova_hyph$residuals,
     main = "Intraradical Hyphae Residuals",
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Q-Q plot
qqnorm(anova_hyph$residuals, main = "Normal Q-Q Plot - Intraradical Hyphae Residuals")
qqline(anova_hyph$residuals, col = "lightblue", lwd = 2)


# ****************************************************
# Filtered dataset ----
# ****************************************************
colonization_filtered <- colonization_2025 %>%
  filter(Treatment_category %in% c(
    "BB_AMF",
    "BB_Bac",
    "BioBead",
    "Fertilizer",
    "Microbes_Liquid",
    "HardyGro",
    "Hydrogel",
    "Neg_ctrl"
  ))

# ****************************************************
# Stacked bar plot ----
# ****************************************************
summary_stacked <- colonization_filtered %>%
  group_by(Treatment_category) %>%
  summarise(
    mean_VS    = mean(Percent_colonization_VS,     na.rm = TRUE),
    mean_Arb   = mean(Percent_colonization_Arb,    na.rm = TRUE),
    mean_Hyph  = mean(Percent_colonization_hyphae, na.rm = TRUE),
    se_total   = sd(Percent_colonization_total, na.rm = TRUE) / sqrt(sum(!is.na(Percent_colonization_total))),
    total_mean = mean(Percent_colonization_total, na.rm = TRUE),
    .groups = "drop"
  )

stacked_long <- summary_stacked %>%
  tidyr::pivot_longer(
    cols      = c(mean_VS, mean_Arb, mean_Hyph),
    names_to  = "Structure",
    values_to = "Mean_Pct"
  ) %>%
  mutate(
    Structure = recode(Structure,
                       "mean_VS"   = "Vesicle/Spore",
                       "mean_Arb"  = "Arbuscule",
                       "mean_Hyph" = "Intraradical Hyphae"
    ),
    Structure = factor(Structure,
                       levels = c("Intraradical Hyphae", "Arbuscule", "Vesicle/Spore")
    )
  )

structure_colors <- c(
  "Intraradical Hyphae" = "#80B1D3",
  "Arbuscule"           = "#BEBADA",
  "Vesicle/Spore"       = "#FB8072"
)

ggplot(stacked_long, aes(x = Treatment_category, y = Mean_Pct, fill = Structure)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(
    data = summary_stacked,
    aes(x = Treatment_category, ymin = total_mean - se_total, ymax = total_mean + se_total),
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

# ****************************************************
# Pearson's correlation ----
# ****************************************************
library(Hmisc)
library(corrplot)

## correlation matrix
master <- read.csv("2025greenhousemaster_8tx.csv")

cor_vars <- master[, c("AMF_colonization", "Arb_colonization", "Hyph_colonization",
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

# ****************************************************
# PCA: AMF to plant responses ----
# ****************************************************
library(vegan)
library(dplyr)
library(tidyr)
library(ggbiplot)
library(ggplot2)

## data standardization
treatments <- master$Treatment_Category

properties <- master[, c("AMF_colonization", "Arb_colonization", "Hyph_colonization",
                         "GrainWeight_g", "Root_biomass_g", "Total_biomass")]

properties_rel <- decostand(properties, method = "max")

##graphic prep
treatment_colors <- c(
  "BioBead_AMF"    = "#8DD3C7",
  "BioBead_Bac"    = "#BC80BD",
  "BioBead_Mixed"        = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Liquid_Microbes"= "#FFFFB3",  
  "HardyGro"       = "#FDB462",
  "Empty_Bead"       = "#80B1D3",
  "Negative"       = "#D9D9D9"    
)

treatment_colors <- c(
  "BioBead_bac"         = "#8DD3C7",
  "BioBead_AMF"         = "#BEBADA",
  "BioBead_Mixed"       = "#FB8072",
  "SL_BioBead_Bac"      = "#FDB462",
  "SL_BioBead_AMF"      = "#80B1D3",
  "SL_BioBead_Mixed"    = "#FFFFB3",
  "Liquid_Mixed"        = "#D9D9D9",
  "HardyGro_Competitor" = "#BC80BD",
  "Empty_Bead"          = "#CCEBC5",
  "SL_Empty_Bead"       = "#FFED6F",
  "Fertilizer"          = "#FC8D62",
  "Negative"            = "#E5C494"
)

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
  "BioBead_AMF"         = "AMF Beads",
  "BioBead_Mixed"       = "AMF Beads",
  "BioBead_Bac"         = "Non-AMF Beads",
  "Empty_Bead"          = "Non-AMF Beads",
  "Liquid_Microbes"        = "Non-Bead Controls",
  "HardyGro" = "Non-Bead Controls",
  "Fertilizer"          = "Non-Bead Controls",
  "Negative"            = "Non-Bead Controls"
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
