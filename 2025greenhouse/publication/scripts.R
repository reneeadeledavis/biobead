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
library(dunn.test)
library(multcomp)
library(FSA)
library(rcompanion)

# ****************************************************
# Load data ----
# ****************************************************
data <- read.csv("publication/2025greenhousemaster.csv")
data$Treatment <- as.factor(data$Treatment) #needed for multcomp package

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
    
    #N2O
    mean_n2o = mean(Fm_Flux_ug.N2O.N_sqmt_hr, na.rm = TRUE),
    sd_n2o   = sd(Fm_Flux_ug.N2O.N_sqmt_hr,   na.rm = TRUE),
    se_n2o   = sd(Fm_Flux_ug.N2O.N_sqmt_hr,   na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

# Display summary statistics
cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)

# ****************************************************
# Plot preparation ----
# ****************************************************

# Set color palette
treatment_labels <- c(
  "BB_AMF"          = "BioBead (AMF)",
  "BB_Bac"          = "BioBead (Bacteria)",
  "BioBead"         = "BioBead (Mix)",
  "Fertilizer"      = "Synthetic Fertilizer",
  "HardyGro"        = "HardyGro",
  "Hydrogel"        = "Empty Bead",
  "Microbes_Liquid" = "Liquid Consortia",
  "Neg_ctrl"        = "No Treatment"
)

treatment_colors <- c(
  # BioBead formulations — blues/teals
  "BB_Bac"          = "#4E9BB5",   # medium blue
  "BB_AMF"          = "#2E7D9C",   # deeper blue
  "BioBead"         = "#1B5E75",   # darkest blue (full mix)
  
  # Carrier/delivery controls — greens
  "Hydrogel"        = "#74C69D",   # empty bead (carrier only)
  "Microbes_Liquid" = "#40916C",   # liquid (same microbes, no bead)
  
  # Commercial comparators — oranges
  "HardyGro"        = "#F4A261",   # commercial biofertilizer
  
  # Conventional/experimental controls — grays & red
  "Fertilizer"      = "#E63946",   # positive control (red stands out)
  "Neg_ctrl"        = "#CCCCCC"    # negative control (neutral gray)
)


# Load Tukey letter function
get_tukey_letters <- function(aov_object) {
  cld_result <- multcomp::cld(
    multcomp::glht(aov_object, linfct = multcomp::mcp(Treatment = "Tukey")),
    decreasing = TRUE   # 'a' = highest group
  )
  # Pull the letter vector into a data frame
  letters_df <- data.frame(
    Treatment = names(cld_result$mcletters$Letters),
    tukey_label = cld_result$mcletters$Letters,
    row.names = NULL
  )
  return(letters_df)
}

get_dunn_letters <- function(data, response_var, group_var = "Treatment", p_adjust = "bh") {
  
  dunn_result <- FSA::dunnTest(
    x      = data[[response_var]],
    g      = data[[group_var]],
    method = p_adjust
  )
  
  dunn_table <- dunn_result$res
  
  letters_out <- rcompanion::cldList(
    formula   = P.adj ~ Comparison,
    data      = dunn_table,
    threshold = 0.05
  )
  
  return(
    letters_out %>%
      rename(Treatment = Group, tukey_label = Letter) %>%
      select(Treatment, tukey_label)
  )
}

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
letters_grainweight <- get_tukey_letters(anova_grainweight)

plot_data_grain <- summary_data %>%
  left_join(letters_grainweight, by = "Treatment")

ggplot(plot_data_grain, aes(x = reorder(Treatment, mean_grain_weight),
                            y = mean_grain_weight,
                            fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_grain_weight - se_grain_weight,
                    ymax = mean_grain_weight + se_grain_weight),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_grain_weight + se_grain_weight + 0.05,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
  scale_x_discrete(labels = treatment_labels) +
  labs(
    title = expression(bold("Grain Weight")),
    y     = expression(bold("Grain Weight (g)")),
    x     = NULL
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y       = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title         = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
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
letters_harvestindex <- get_tukey_letters(anova_harvestindex)

plot_data_harvestindex <- summary_data %>%
  left_join(letters_harvestindex, by = "Treatment")

ggplot(plot_data_harvestindex, aes(x = reorder(Treatment, mean_harvest_index), 
                         y = mean_harvest_index,
                         fill = Treatment)) +  
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +  
  scale_fill_manual(values = treatment_colors) +  
  geom_errorbar(aes(ymin = mean_harvest_index - se_harvest_index, 
                    ymax = mean_harvest_index + se_harvest_index),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_harvest_index + se_harvest_index + 0.05,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
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
### Kruskal-Wallis ----

# NH4
kruskal.test(NH4.N_ugml ~ Treatment, data = data)
# Kruskal-Wallis chi-squared = 13.556, df = 7, p-value = 0.05966
dunn.test(data$NH4.N_ugml, data$Treatment, method = "bh")

# NO3
kruskal.test(NO3.N_ugml ~ Treatment, data = data)
# Kruskal-Wallis chi-squared = 16.553, df = 7, p-value = 0.02052
dunn.test(data$NO3.N_ugml, data$Treatment, method = "bh")

# PO4
kruskal.test(PO4.P_ugml ~ Treatment, data = data)
# Kruskal-Wallis chi-squared = 13.503, df = 7, p-value = 0.06076
dunn.test(data$PO4.P_ugml, data$Treatment, method = "bh")

### Dunn CLD letters ----
letters_nh4 <- get_dunn_letters(data, "NH4.N_ugml") %>% mutate(Nutrient = "NH4-N")
letters_no3 <- get_dunn_letters(data, "NO3.N_ugml") %>% mutate(Nutrient = "NO3-N")
letters_po4 <- get_dunn_letters(data, "PO4.P_ugml") %>% mutate(Nutrient = "PO4-P")

leachate_letters <- bind_rows(letters_nh4, letters_no3, letters_po4)

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
  ) %>%
  left_join(leachate_letters, by = c("Treatment", "Nutrient"))

ggplot(leachate_long, aes(x = Treatment, y = Mean_ugml, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean_ugml - SE, ymax = Mean_ugml + SE),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = Mean_ugml + SE, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Nutrient, scales = "free_y") +
  scale_fill_manual(values = treatment_colors) +
  scale_x_discrete(labels = treatment_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
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

## N2O ----
### Kruskal-Wallis ----
kruskal.test(Fm_Flux_ug.N2O.N_sqmt_hr ~ Treatment, data = data)
# Kruskal-Wallis chi-squared = 14.857, df = 7, p-value = 0.03787
dunn.test(data$Fm_Flux_ug.N2O.N_sqmt_hr, data$Treatment, method = "bh")

### Bar plot ----
ggplot(summary_data, aes(x = reorder(Treatment, mean_n2o),
                        y = mean_n2o,
                        fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_n2o - se_n2o,
                    ymax = mean_n2o + se_n2o),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("N"[2]*"O-N Flux by Treatment")),
    y     = expression(bold("µg N"[2]*"O-N / m² / hr (Mean ± SE)")),
    x     = NULL
  ) +
  theme_bw(base_size = 12) +
  guides(fill = "none") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x       = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y       = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title         = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
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
letters_colonization <- get_tukey_letters(anova_colonization)

plot_data_colonization <- summary_data %>%
  left_join(letters_colonization, by = "Treatment")

ggplot(plot_data_colonization, aes(x = reorder(Treatment, mean_pct_total),
                         y = mean_pct_total,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_total - se_pct_total,
                    ymax = mean_pct_total + se_pct_total),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_pct_total + se_pct_total + 2,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) + #raises y to fit letters
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
letters_arbuscule <- get_tukey_letters(anova_arbuscule)

plot_data_arbuscule <- summary_data %>%
  left_join(letters_arbuscule, by = "Treatment")

ggplot(plot_data_arbuscule, aes(x = reorder(Treatment, mean_pct_arb),
                         y = mean_pct_arb,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_arb - se_pct_arb,
                    ymax = mean_pct_arb + se_pct_arb),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_pct_arb + se_pct_arb + 0.5,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
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
letters_hyph <- get_tukey_letters(anova_hyph)

plot_data_hyph <- summary_data %>%
  left_join(letters_hyph, by = "Treatment")

ggplot(plot_data_hyph, aes(x = reorder(Treatment, mean_pct_hyph),
                         y = mean_pct_hyph,
                         fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_pct_hyph - se_pct_hyph,
                    ymax = mean_pct_hyph + se_pct_hyph),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_pct_hyph + se_pct_hyph + 1,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
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
