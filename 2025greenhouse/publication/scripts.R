# ****************************************************
# Load packages ----
# ****************************************************
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)
library(vegan)
library(ggbiplot)
library(multcomp)
library(FSA)
library(rcompanion)
library(tidytext)

# ****************************************************
# Load data ----
# ****************************************************
data <- read.csv("publication/2025greenhousemaster.csv", fileEncoding = "UTF-8-BOM")
data$Treatment <- as.factor(data$Treatment) # needed for multcomp package

tissue <- read.csv("publication/OSU.csv", fileEncoding = "UTF-8-BOM")
tissue$Treatment <- as.factor(tissue$Treatment)
tissue$Tissue    <- as.factor(tissue$Tissue) # making everything a factor for multcomp

# create new N (g) and P (g) columns from %N and %P
tissue <- tissue %>%
  mutate(
    P_g = Biomass_g * (P / 100), # original data is %, not proportion
    N_g = Biomass_g * (N / 100)
  )

# ****************************************************
# Summary statistics ----
# ****************************************************
## Whole plant data ----
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
    
    
    # AMF colonization
    mean_pct_total = mean(Hyph_colonization + Arb_colonization, na.rm = TRUE),
    sd_pct_total   = sd(Hyph_colonization + Arb_colonization, na.rm = TRUE),
    se_pct_total   = sd(Hyph_colonization + Arb_colonization, na.rm = TRUE) / sqrt(n()),
    
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

## Tissue nutrients ----
# need to group by Tissue AND Treatment so each facet gets its own mean/SE
summary_tissue <- tissue %>%
  group_by(Tissue, Treatment) %>%
  summarise(
    mean_C = mean(C, na.rm = TRUE),
    se_C   = sd(C, na.rm = TRUE) / sqrt(n()),
    
    mean_N = mean(N, na.rm = TRUE),
    se_N   = sd(N, na.rm = TRUE) / sqrt(n()),
    
    mean_CN = mean(CN_ratio, na.rm = TRUE),
    se_CN   = sd(CN_ratio, na.rm = TRUE) / sqrt(n()),
    
    mean_prot = mean(Crude_protein, na.rm = TRUE),
    se_prot   = sd(Crude_protein, na.rm = TRUE) / sqrt(n()),
    
    mean_P = mean(P, na.rm = TRUE),
    se_P   = sd(P, na.rm = TRUE) / sqrt(n()),
    
    mean_P_g = mean(P_g, na.rm = TRUE),
    se_P_g   = sd(P_g, na.rm = TRUE) / sqrt(n()),
    
    mean_N_g = mean(N_g, na.rm = TRUE),
    se_N_g   = sd(N_g, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = "drop"
  )

# create new N (g) and P (g) columns from %N and %P
tissue <- tissue %>%
  mutate(
    P_g = Biomass_g * (P / 100), #original data is %, not proportion
    N_g = Biomass_g * (N / 100)
  )

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


# Helper functions ----
get_tukey_letters <- function(aov_object) {
  cld_result <- multcomp::cld(
    multcomp::glht(aov_object, linfct = multcomp::mcp(Treatment = "Tukey")),
    decreasing = TRUE, #a = highest group
    Letters    = letters, # forces b and c, otherwise they might be dropped
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
      dplyr::rename(Treatment = Group, tukey_label = Letter) %>%
      dplyr::select(Treatment, tukey_label)        # explicit namespace
  )
}

# Combine per-tissue Tukey letters into one table for faceted plots
combine_tissue_letters <- function(root, abg, grain) {
  bind_rows(
    root  %>% mutate(Tissue = "Root"),
    abg   %>% mutate(Tissue = "ABG"),
    grain %>% mutate(Tissue = "Grain")
  ) %>%
    mutate(
      Treatment   = as.character(Treatment), # ensure character for join
      tukey_label = trimws(tukey_label)      # strip cld() padding spaces
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
tukey_grainweight

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

### Tukey letters ----
letters_grainweight <- get_tukey_letters(anova_grainweight)

### Bar plot ----
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

## Harvest index ----
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

### Tukey letters ----
letters_harvestindex <- get_tukey_letters(anova_harvestindex)

### Bar plot ----
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
  geom_text(aes(y = mean_harvest_index + se_harvest_index + 0.015,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
  scale_x_discrete(labels = treatment_labels) +
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

# ****************************************************
# Part 2: Plant tissue nutrients ----
# ****************************************************
## Crude protein ----
# Crude protein in roots
anova_prot_r <- aov(Crude_protein ~ Treatment, data = r)
summary(anova_prot_r)

tukey_prot_r <- TukeyHSD(anova_prot_r)
tukey_prot_r

shapiro.test(anova_prot_r$residuals)
# W = 0.97314, p-value = 0.7445

# Crude protein in above ground biomass
anova_prot_abg <- aov(Crude_protein ~ Treatment, data = abg)
summary(anova_prot_abg)

tukey_prot_abg <- TukeyHSD(anova_prot_abg)
tukey_prot_abg

shapiro.test(anova_prot_abg$residuals)
# W = 0.8867, p-value = 0.01137

# Crude protein in grain
anova_prot_g <- aov(Crude_protein ~ Treatment, data = g)
summary(anova_prot_g)

tukey_prot_g <- TukeyHSD(anova_prot_g)
tukey_prot_g

shapiro.test(anova_prot_g$residuals)
# W = 0.97231, p-value = 0.7242

### Tukey letters ----
letters_prot_r   <- get_tukey_letters(anova_prot_r)
letters_prot_abg <- get_tukey_letters(anova_prot_abg)
letters_prot_g   <- get_tukey_letters(anova_prot_g)
letters_prot     <- combine_tissue_letters(letters_prot_r, letters_prot_abg, letters_prot_g)

### Faceted bar plot (all tissues) ----
tissue_prot <- summary_tissue %>%
  left_join(letters_prot, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_prot, Tissue)
  )

ggplot(tissue_prot, aes(x = Treatment_ordered, y = mean_prot, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_prot - se_prot, ymax = mean_prot + se_prot),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_prot + se_prot, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Crude Protein Estimation by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean Crude Protein Estimation (%)")),
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

### Grain-only bar plot ----
plot_data_prot_g <- summary_tissue %>%
  filter(Tissue == "Grain") %>%
  left_join(letters_prot_g, by = "Treatment")

ggplot(plot_data_prot_g, aes(x = reorder(Treatment, mean_prot), 
                             y = mean_prot,
                             fill = Treatment)) +  
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +  
  scale_fill_manual(values = treatment_colors) +  
  geom_errorbar(aes(ymin = mean_prot - se_prot, 
                    ymax = mean_prot + se_prot),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_prot + se_prot + 0.3,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
  scale_x_discrete(labels = treatment_labels) +
  labs(
    title = expression(bold("Crude Protein Estimation in Grain")),
    y = expression(bold("Crude Protein Estimation (%)")),
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

## C conc (%C) ----
### ANOVA ----
# C in roots
anova_C_r <- aov(C ~ Treatment, data = r)
summary(anova_C_r)

tukey_C_r <- TukeyHSD(anova_C_r)
tukey_C_r

shapiro.test(anova_C_r$residuals)
# W = 0.76817, p-value = 9.286e-05

# C in above ground biomass
anova_C_abg <- aov(C ~ Treatment, data = abg)
summary(anova_C_abg)

tukey_C_abg <- TukeyHSD(anova_C_abg)
tukey_C_abg

shapiro.test(anova_C_abg$residuals)
# W = 0.9834, p-value = 0.9494

# C in grain
anova_C_g <- aov(C ~ Treatment, data = g)
summary(anova_C_g)

tukey_C_g <- TukeyHSD(anova_C_g)
tukey_C_g

shapiro.test(anova_C_g$residuals)
# W = 0.97803, p-value = 0.857

### Tukey letters ----
letters_C <- combine_tissue_letters(get_tukey_letters(anova_C_r),
                                    get_tukey_letters(anova_C_abg),
                                    get_tukey_letters(anova_C_g))

### Faceted bar plot ----
tissue_C <- summary_tissue %>%
  left_join(letters_C, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_C, Tissue)
  )

ggplot(tissue_C, aes(x = Treatment_ordered, y = mean_C, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_C - se_C, ymax = mean_C + se_C),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_C + se_C, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Carbon Concentration by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean C (%)")),
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

## N conc (%N) ----
### ANOVA ----
# N in roots
anova_N_r <- aov(N ~ Treatment, data = r)
summary(anova_N_r)

tukey_N_r <- TukeyHSD(anova_N_r)
tukey_N_r

shapiro.test(anova_N_r$residuals)
# W = 0.97314, p-value = 0.7445

# N in above ground biomass
anova_N_abg <- aov(N ~ Treatment, data = abg)
summary(anova_N_abg)

tukey_N_abg <- TukeyHSD(anova_N_abg)
tukey_N_abg

shapiro.test(anova_N_abg$residuals)
# W = 0.8867, p-value = 0.01137 - not great 

# Quantile-Quantile plots
qqnorm(anova_N_abg$residuals, main = "Normal Q-Q Plot - %N in Above Ground Biomass Residuals")
qqline(anova_N_abg$residuals, col = "lightblue", lwd = 2)

kruskal.test(N ~ Treatment, data = abg)

# N in grain
anova_N_g <- aov(N ~ Treatment, data = g)
summary(anova_N_g)

tukey_N_g <- TukeyHSD(anova_N_g)
tukey_N_g

shapiro.test(anova_N_g$residuals)
# W = 0.97231, p-value = 0.7242

### Tukey letters ----
letters_N <- combine_tissue_letters(get_tukey_letters(anova_N_r),
                                    get_tukey_letters(anova_N_abg),
                                    get_tukey_letters(anova_N_g))

### Faceted bar plot ----
tissue_N <- summary_tissue %>%
  left_join(letters_N, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_N, Tissue)
  )

ggplot(tissue_N, aes(x = Treatment_ordered, y = mean_N, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_N - se_N, ymax = mean_N + se_N),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_N + se_N, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Nitrogen Concentration by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean N (%)")),
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

## N mass (g) ----
### ANOVA ----
# N in roots
anova_N_g_r <- aov(N_g ~ Treatment, data = r)
summary(anova_N_g_r)

tukey_N_g_r <- TukeyHSD(anova_N_g_r)
tukey_N_g_r 

shapiro.test(anova_N_g_r$residuals)
# W = 0.98, p-value = 0.8958

# N in above ground biomass
anova_N_g_abg <- aov(N_g ~ Treatment, data = abg)
summary(anova_N_g_abg)

tukey_N_g_abg <- TukeyHSD(anova_N_g_abg)
tukey_N_g_abg

shapiro.test(anova_N_g_abg$residuals)
# W = 0.87699, p-value = 0.007231; data not normal

# Quantile-Quantile plots
qqnorm(anova_N_g_abg$residuals, main = "Normal Q-Q Plot - N mass in Above Ground Biomass Residuals")
qqline(anova_N_g_abg$residuals, col = "lightblue", lwd = 2)

kruskal.test(N_g ~ Treatment, data = abg)

# N in grain
anova_N_g_g <- aov(N_g ~ Treatment, data = g)
summary(anova_N_g_g)

tukey_N_g_g <- TukeyHSD(anova_N_g_g)
tukey_N_g_g

shapiro.test(anova_N_g_g$residuals)
# W = 0.98306, p-value = 0.9449

### Tukey letters ----
letters_N_g <- combine_tissue_letters(get_tukey_letters(anova_N_g_r),
                                      get_tukey_letters(anova_N_g_abg),
                                      get_tukey_letters(anova_N_g_g))

### Faceted bar plot ----
tissue_N_g <- summary_tissue %>%
  left_join(letters_N_g, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_N_g, Tissue)
  )

ggplot(tissue_N_g, aes(x = Treatment_ordered, y = mean_N_g, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_N_g - se_N_g, ymax = mean_N_g + se_N_g),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_N_g + se_N_g, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Total N mass by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean N mass (g)")),
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

## C:N ratio ----
### ANOVA ----
# C:N in roots
anova_CN_r <- aov(CN_ratio ~ Treatment, data = r)
summary(anova_CN_r)

tukey_CN_r <- TukeyHSD(anova_CN_r)
tukey_CN_r

shapiro.test(anova_CN_r$residuals)
# W = 0.94115, p-value = 0.173

# C:N in above ground biomass
anova_CN_abg <- aov(CN_ratio ~ Treatment, data = abg)
summary(anova_CN_abg)

tukey_CN_abg <- TukeyHSD(anova_CN_abg)
tukey_CN_abg

shapiro.test(anova_CN_abg$residuals)
# W = 0.97644, p-value = 0.8223

# C:N in grain
anova_CN_g <- aov(CN_ratio ~ Treatment, data = g)
summary(anova_CN_g)

tukey_CN_g <- TukeyHSD(anova_CN_g)
tukey_CN_g

shapiro.test(anova_CN_g$residuals)
# W = 0.95896, p-value = 0.4178

### Tukey letters ----
letters_CN <- combine_tissue_letters(get_tukey_letters(anova_CN_r),
                                     get_tukey_letters(anova_CN_abg),
                                     get_tukey_letters(anova_CN_g))

### Faceted bar plot ----
tissue_CN <- summary_tissue %>%
  left_join(letters_CN, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_CN, Tissue)
  )

ggplot(tissue_CN, aes(x = Treatment_ordered, y = mean_CN, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_CN - se_CN, ymax = mean_CN + se_CN),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_CN + se_CN, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("C:N Ratio by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean C:N Ratio")),
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

## P conc (%P) ----
### ANOVA ----
# P in roots
anova_P_r <- aov(P ~ Treatment, data = r)
summary(anova_P_r)

tukey_P_r <- TukeyHSD(anova_P_r)
tukey_P_r

shapiro.test(anova_P_r$residuals)
# W = 0.94091, p-value = 0.171

# P in above ground biomass
anova_P_abg <- aov(P ~ Treatment, data = abg)
summary(anova_P_abg)

tukey_P_abg <- TukeyHSD(anova_P_abg)
tukey_P_abg

shapiro.test(anova_P_abg$residuals)
# W = 0.94206, p-value = 0.1813

# P in grain
anova_P_g <- aov(P ~ Treatment, data = g)
summary(anova_P_g)

tukey_P_g <- TukeyHSD(anova_P_g)
tukey_P_g

shapiro.test(anova_P_g$residuals)
# W = 0.96745, p-value = 0.6045

### Tukey letters ----
letters_P <- combine_tissue_letters(get_tukey_letters(anova_P_r),
                                    get_tukey_letters(anova_P_abg),
                                    get_tukey_letters(anova_P_g))

### Faceted bar plot ----
tissue_P <- summary_tissue %>%
  left_join(letters_P, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_P, Tissue)
  )

ggplot(tissue_P, aes(x = Treatment_ordered, y = mean_P, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_P - se_P, ymax = mean_P + se_P),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_P + se_P, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Phosphorus Concentration by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean P (%)")),
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

## P mass (g) ----
### ANOVA ----
# P in roots
anova_P_g_r <- aov(P_g ~ Treatment, data = r)
summary(anova_P_g_r)

tukey_P_g_r <- TukeyHSD(anova_P_g_r)
tukey_P_g_r 

shapiro.test(anova_P_g_r$residuals)
# W = 0.96397, p-value = 0.5233

# P in above ground biomass
anova_P_g_abg <- aov(P_g ~ Treatment, data = abg)
summary(anova_P_g_abg)

tukey_P_g_abg <- TukeyHSD(anova_P_g_abg)
tukey_P_g_abg

shapiro.test(anova_P_g_abg$residuals)
# W = 0.96034, p-value = 0.4452

# P in grain
anova_P_g_g <- aov(P_g ~ Treatment, data = g)
summary(anova_P_g_g)

tukey_P_g_g <- TukeyHSD(anova_P_g_g)
tukey_P_g_g

shapiro.test(anova_P_g_g$residuals)
# W = 0.94781, p-value = 0.2427

### Tukey letters ----
letters_P_g <- combine_tissue_letters(get_tukey_letters(anova_P_g_r),
                                      get_tukey_letters(anova_P_g_abg),
                                      get_tukey_letters(anova_P_g_g))

### Faceted bar plot ----
tissue_P_g <- summary_tissue %>%
  left_join(letters_P_g, by = c("Tissue", "Treatment")) %>%
  mutate(
    Treatment         = factor(Treatment, levels = names(treatment_colors)),
    Treatment_display = recode(as.character(Treatment), !!!treatment_labels),
    Treatment_ordered = reorder_within(Treatment_display, mean_P_g, Tissue)
  )

ggplot(tissue_P_g, aes(x = Treatment_ordered, y = mean_P_g, fill = Treatment)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_P_g - se_P_g, ymax = mean_P_g + se_P_g),
                width = 0.2, linewidth = 0.6, color = "grey30") +
  geom_text(aes(y = mean_P_g + se_P_g, label = tukey_label),
            size = 3.5, vjust = -0.4) +
  facet_wrap(~ Tissue, scales = "free", labeller = as_labeller(tissue_labels)) +
  scale_x_reordered() +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = expression(bold("Total P mass by Plant Tissue")),
    x     = NULL,
    y     = expression(bold("Mean P mass (g)")),
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

## Whole plant P ----
### Data prep ----
whole_plant_p <- tissue %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(Whole_Plant_P_g = sum(P_g), .groups = "drop")

summary_whole_plant_p <- whole_plant_p %>%
  group_by(Treatment) %>%
  summarise(
    mean_whole_plant_P = mean(Whole_Plant_P_g, na.rm = TRUE),
    se_whole_plant_P   = sd(Whole_Plant_P_g, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

### ANOVA ----
anova_whole_plant_P <- aov(Whole_Plant_P_g ~ Treatment, data = whole_plant_p)
summary(anova_whole_plant_P)

tukey_whole_plant_P <- TukeyHSD(anova_whole_plant_P)
tukey_whole_plant_P

# Shapiro Wilk test for normality
shapiro.test(anova_whole_plant_P$residuals)
# W = 0.97366, p-value = 0.757 - data normally distributed

# Residuals distribution
hist(anova_whole_plant_P$residuals, 
     main = "Whole Plant P Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_whole_plant_P$residuals, main = "Normal Q-Q Plot - Whole Plant P Residuals")
qqline(anova_whole_plant_P$residuals, col = "lightblue", lwd = 2)

### Tukey letters ----
# treatment-level (not tissue-level), so no combine_tissue_letters()
letters_whole_plant_P <- get_tukey_letters(anova_whole_plant_P)

### Bar plot ----
plot_data_whole_plant_P <- summary_whole_plant_p %>%
  left_join(letters_whole_plant_P, by = "Treatment")

ggplot(plot_data_whole_plant_P, aes(x = reorder(Treatment, mean_whole_plant_P),
                                    y = mean_whole_plant_P,
                                    fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_whole_plant_P - se_whole_plant_P,
                    ymax = mean_whole_plant_P + se_whole_plant_P),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_whole_plant_P + se_whole_plant_P + 0.001,  # nudge above error bar
                label = tukey_label),
            size = 4, vjust = 0) +
  scale_x_discrete(labels = treatment_labels) +
  labs(
    title = expression(bold("Whole Plant P")),
    y     = expression(bold("Phosphorous (g)")),
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

## P harvest index ----
# grain P (g) / whole-plant P (g)
### Data prep ----
p_harvest <- tissue %>%
  group_by(Sample_ID, Treatment) %>%
  summarise(
    Grain_P_g       = sum(P_g[Tissue == "Grain"]),
    Whole_Plant_P_g = sum(P_g),
    .groups = "drop"
  ) %>%
  mutate(P_Harvest_Index = Grain_P_g / Whole_Plant_P_g)

summary_p_harvest <- p_harvest %>%
  group_by(Treatment) %>%
  summarise(
    mean_P_HI = mean(P_Harvest_Index, na.rm = TRUE),
    sd_P_HI   = sd(P_Harvest_Index, na.rm = TRUE),
    se_P_HI   = sd(P_Harvest_Index, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )
print(summary_p_harvest)

### ANOVA ----
anova_p_harvest <- aov(P_Harvest_Index ~ Treatment, data = p_harvest)
summary(anova_p_harvest)

tukey_p_harvest <- TukeyHSD(anova_p_harvest)
tukey_p_harvest

# Shapiro-Wilk
shapiro.test(anova_p_harvest$residuals)
# W = 0.97274, p-value = 0.7345

# Residual diagnostics
hist(anova_p_harvest$residuals,
     main = "P Harvest Index Residuals",
     xlab = "Residuals", col = "lightblue", breaks = 20)
qqnorm(anova_p_harvest$residuals, main = "Normal Q-Q Plot - P Harvest Index Residuals")
qqline(anova_p_harvest$residuals, col = "lightblue", lwd = 2)

### Tukey letters ----
letters_p_harvest <- get_tukey_letters(anova_p_harvest)

### Bar plot ----
plot_data_p_harvest <- summary_p_harvest %>%
  left_join(letters_p_harvest, by = "Treatment")

ggplot(plot_data_p_harvest, aes(x = reorder(Treatment, mean_P_HI),
                                y = mean_P_HI,
                                fill = Treatment)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = treatment_colors) +
  geom_errorbar(aes(ymin = mean_P_HI - se_P_HI,
                    ymax = mean_P_HI + se_P_HI),
                width = 0.3, linewidth = 0.5) +
  geom_text(aes(y = mean_P_HI + se_P_HI + 0.02, label = tukey_label),
            size = 4, vjust = 0) +
  scale_x_discrete(labels = treatment_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = expression(bold("Phosphorus Harvest Index")),
    y     = expression(bold("P Harvest Index (grain P / total plant P)")),
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

# ****************************************************
# Part 3: AMF establishment ----
# ****************************************************

## Stacked bar plot ----
summary_stacked <- data %>%
  group_by(Treatment) %>%
  summarise(
    mean_Arb   = mean(Arb_colonization,  na.rm = TRUE),
    mean_Hyph  = mean(Hyph_colonization, na.rm = TRUE),
    se_total   = sd(Hyph_colonization + Arb_colonization, na.rm = TRUE) / sqrt(sum(!is.na(Hyph_colonization))),
    total_mean = mean(Hyph_colonization + Arb_colonization, na.rm = TRUE),
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

ggplot(stacked_long, aes(x = reorder(Treatment, Mean_Pct, FUN = sum), y = Mean_Pct, fill = Structure)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_errorbar(
    data = summary_stacked,
    aes(x = reorder(Treatment, total_mean), ymin = total_mean - se_total, ymax = total_mean + se_total),    inherit.aes = FALSE,
    width = 0.2, linewidth = 0.6, color = "grey30"
  ) +
  scale_fill_manual(values = structure_colors) +
  scale_x_discrete(labels = treatment_labels) +
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
  scale_x_discrete(labels = treatment_labels) +
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
  scale_x_discrete(labels = treatment_labels) +
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
  scale_x_discrete(labels = treatment_labels) +
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
# Part 4: Correlations ----
# ****************************************************

## Data prep ----
## Pool composite ----
# G{trt}-{k} = master plants 2k-1 and 2k (G1-1 = T1-1 + T1-2); plant 7 was not analyzed
data_pool <- data %>%
  tidyr::separate(Sample_ID, into = c("trt_num", "plant"), sep = "-", remove = FALSE) %>%
  mutate(trt_num = sub("^T", "", trt_num),
         plant   = as.integer(plant)) %>%
  filter(plant <= 6) %>%
  mutate(Pool_ID = paste0("G", trt_num, "-", ceiling(plant / 2))) %>% #this is how the pool is created
  group_by(Pool_ID, Treatment) %>%
  summarise(
    # order is important - HI first: the next lines overwrite columns of the same name inside summarise()
    Harvest_Index    = sum(GrainWeight_g) / sum(Total_biomass - Root_biomass_g),
    AMF_colonization = mean(AMF_colonization),
    Total_biomass    = mean(Total_biomass),
    n_plants         = n(),
    .groups = "drop"
  )

tissue$Pool_ID <- sub("^T", "G", tissue$Sample_ID)  # matches with OSU data, creates a new column from the result

tissue_wide <- tissue %>%
  dplyr::select(Pool_ID, Tissue, N, P, P_g) %>%
  pivot_wider(names_from = Tissue, values_from = c(N, P, P_g), names_sep = "_") %>%
  mutate(
    Whole_Plant_P_g = (P_g_Root + P_g_ABG + P_g_Grain) / 2,    # per plant; each pool = 2 plants
    P_Harvest_Index = P_g_Grain / (P_g_Root + P_g_ABG + P_g_Grain)
  )

composite <- inner_join(data_pool, tissue_wide, by = "Pool_ID")
stopifnot(nrow(composite) == 24, all(composite$n_plants == 2))

## Pearson's correlation ----
cor_vars <- composite[, c("N_Root", "N_ABG", "N_Grain",          # %N by tissue
                          "P_Root", "P_ABG", "P_Grain",          # %P by tissue
                          "Whole_Plant_P_g", "P_Harvest_Index",
                          "AMF_colonization",                    # total AMF
                          "Total_biomass", "Harvest_Index")]

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

## Correlation plot ----
# from Korena Mafune, unedited/pasted in
# 1. Compute the correlation matrix
cor.n.ss <- Lettuce[,9:15]
M <- cor(cor.n.ss)

# 2. Create the basic plot
corrplot(M, method = "circle") 
corrplot.mixed(M, order = 'AOE')
testRes = cor.mtest(M, conf.level = 0.95)

## specialized the insignificant value according to the significant level
corrplot(M, p.mat = testRes$p, sig.level = 0.05)

ggscatter(dataset, x = "Arbuscule", y = "Yield_kg", color = "black", shape = 21, size = 2,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Overall establishment (%)", ylab = "Plot Yield (kg)",
          add.params = list(color = "black", fill = "darkgrey"))

## PCA ----
## data standardization
treatments <- composite$Treatment

properties <- composite[, c("N_Root", "N_ABG", "N_Grain",          # %N by tissue
                         "P_Root", "P_ABG", "P_Grain",          # %P by tissue
                         "Whole_Plant_P_g", "P_Harvest_Index",
                         "AMF_colonization",                    # total AMF
                         "Total_biomass", "Harvest_Index")]

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
  scale_color_manual(values = treatment_colors, 
                     labels = treatment_labels,
                     breaks = names(treatment_colors),
                     name = "Treatment") +
  scale_x_continuous(expand = expansion(mult = 0.5)) +
  scale_y_continuous(expand = expansion(mult = 0.3)) +
  labs(colour = "Treatment",
       x = "PC1 (58.2%)",   # update from summary()
       y = "PC2 (13.3%)",
       caption = "PCA: AMF colonization and plant response variables") +
  theme_bw()

pca_plot
