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
wheatyield_2025 <- read.csv("yield/Wheat_Biomass_2025.csv", header = TRUE) 

# ****************************************************
# Summary statistics ----
# ****************************************************
summary_data <- wheatyield_2025 %>%
  group_by(Treatment) %>%
  summarise(
    
    # GrainHeads 
    mean_grainhead = mean(GrainHeads, na.rm = TRUE),
    sd_grainhead = sd(GrainHeads, na.rm = TRUE),
    se_grainhead = sd(GrainHeads, na.rm = TRUE) / sqrt(n()),
    
    # ABG_straw 
    mean_abg_straw = mean(ABG_straw, na.rm = TRUE),
    sd_abg_straw = sd(ABG_straw, na.rm = TRUE),
    se_abg_straw = sd(ABG_straw, na.rm = TRUE) / sqrt(n()),
    
    # GrainHeadWeight
    mean_grainhead_weight = mean(GrainHeadWeight_g, na.rm = TRUE),
    sd_grainhead_weight = sd(GrainHeadWeight_g, na.rm = TRUE),
    se_grainhead_weight = sd(GrainHeadWeight_g, na.rm = TRUE) / sqrt(n()),
    
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
    
    # Root-Shoot ratio
    mean_root_shoot = mean(Root_Shoot, na.rm = TRUE),
    sd_root_shoot = sd(Root_Shoot, na.rm = TRUE),
    se_root_shoot = sd(Root_Shoot, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

# Display summary statistics
cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)

# Set color palette
treatment_colors <- c(
  "BB_AMF"    = "#8DD3C7",
  "BB_Bac"    = "#BC80BD",
  "BioBead"        = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Microbes_Liquid"= "#FFFFB3",  
  "Hardygro"       = "#FDB462",
  "Hydrogel"       = "#80B1D3",
  "Neg_ctrl"       = "#D9D9D9"    
)

# ****************************************************
# Grain weight bar plot ----
# ****************************************************
## 8 tx ----
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

## 12 tx ----
summary_data <- master %>%
  group_by(Treatment_Category) %>%
  summarise(
    
    # Total_biomass
    mean_Total_biomass = mean(Total_biomass, na.rm = TRUE),
    sd_Total_biomass = sd(Total_biomass, na.rm = TRUE),
    se_Total_biomass = sd(Total_biomass, na.rm = TRUE) / sqrt(n()),
    
    # Grain_weight
    mean_Grain_weight = mean(Grain_weight, na.rm = TRUE),
    sd_Grain_weight = sd(Grain_weight, na.rm = TRUE),
    se_Grain_weight = sd(Grain_weight, na.rm = TRUE) / sqrt(n()),
    
    # Harvest_Index
    mean_harvest_index = mean(Harvest_Index, na.rm = TRUE),
    sd_harvest_index = sd(Harvest_Index, na.rm = TRUE),
    se_harvest_index = sd(Harvest_Index, na.rm = TRUE) / sqrt(n()),
    
    # Root_biomass
    mean_root_biomass = mean(Root_biomass, na.rm = TRUE),
    sd_root_biomass = sd(Root_biomass, na.rm = TRUE),
    se_root_biomass = sd(Root_biomass, na.rm = TRUE) / sqrt(n()),
    
    # Root-Shoot ratio
    mean_root_shoot = mean(Root_Shoot, na.rm = TRUE),
    sd_root_shoot = sd(Root_Shoot, na.rm = TRUE),
    se_root_shoot = sd(Root_Shoot, na.rm = TRUE) / sqrt(n()),
    
    n = n(),
    .groups = 'drop'
  )

# Display summary statistics
cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)

# Set color palette
treatment_colors <- c(
  "BioBead_bac"         = "#8DD3C7",
  "BioBead_AMF"         = "#BEBADA",
  "BioBead_Mixed"       = "#FB8072",
  "SL_BioBead_Mixed"    = "#FFFFB3",
  "SL_BioBead_Bac"      = "#FDB462",
  "SL_BioBead_AMF"      = "#80B1D3",
  "Liquid_Mixed"        = "#E5C494",
  "HardyGro_Competitor" = "#BC80BD",
  "Empty_Bead"          = "#CCEBC5",
  "SL_Empty_Bead"       = "#FFED6F",
  "Fertilizer"          = "#FC8D62",
  "Negative"            = "#D9D9D9"
)

ggplot(summary_data, aes(x = reorder(Treatment_Category, mean_Grain_weight), 
                         y = mean_Grain_weight,
                         fill = Treatment_Category)) +  
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.3) +  
  scale_fill_manual(values = treatment_colors) +  
  geom_errorbar(aes(ymin = mean_Grain_weight - se_Grain_weight, 
                    ymax = mean_Grain_weight + se_Grain_weight),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("Grain Weight")),
    x = NULL,
    y = expression(bold("Grain Weight (g)"))
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
# Harvest index bar plot ----
# ****************************************************

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

# ****************************************************
# Stacked bar plot ----
# ****************************************************

# per-pot totals and husk, then means + SE
stack_summary <- wheatyield_2025 %>%
  mutate(Husk_g  = GrainHeadWeight_g - GrainWeight_g,
         Total_g = Root_biomass_g + ABG_straw + Husk_g + GrainWeight_g) %>%
  group_by(Treatment) %>%
  summarise(
    Root      = mean(Root_biomass_g, na.rm = TRUE),
    Straw     = mean(ABG_straw,      na.rm = TRUE),
    Husk      = mean(Husk_g,         na.rm = TRUE),
    Grain     = mean(GrainWeight_g,  na.rm = TRUE),
    se_total  = sd(Total_g, na.rm = TRUE) / sqrt(sum(!is.na(Total_g))),
    total_mean = mean(Total_g, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: pivot to long for stacking (keep se_total and total_mean as sidebar cols)
stack_data <- stack_summary %>%
  tidyr::pivot_longer(
    cols      = c(Root, Straw, Husk, Grain),
    names_to  = "Component",
    values_to = "Biomass_g"
  ) %>%
  mutate(Component = factor(Component, levels = c("Grain", "Husk", "Straw", "Root")))

# Color palette
component_colors <- c(
  "Root"  = "#8B5E3C",
  "Straw" = "#A8C97F",
  "Husk"  = "#D4B483",
  "Grain" = "#E8834E"
)

# Treatment order
treatment_order <- c("BB_Bac", "BB_AMF", "BioBead", "Fertilizer",
                     "Hardygro", "Hydrogel", "Microbes_Liquid", "Neg_ctrl")

stack_data <- stack_data %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order))

# Step 3: SE bar data (one row per treatment, positioned at total bar top)
se_data <- stack_summary %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order))

# Plot
biomasspartion_plot <- ggplot(stack_data, aes(x = Treatment, y = Biomass_g, fill = Component)) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.3) +
  # SE error bar at top of total bar
  geom_errorbar(
    data = se_data,
    aes(x = Treatment, ymin = total_mean - se_total, ymax = total_mean + se_total),
    inherit.aes = FALSE,
    width = 0.2, linewidth = 0.6, color = "grey30"
  ) +
  scale_fill_manual(
    values = component_colors,
    guide  = guide_legend(reverse = TRUE)                          # flips display to Root on top
  ) +
  scale_x_discrete(labels = c(
    "BB_Bac"         = "BB_Bac",
    "BB_AMF"         = "BB_AMF",
    "BioBead"         = "BioBead",
    "Fertilizer"      = "Fertilizer",
    "Hardygro"        = "Hardygro",
    "Hydrogel"        = "Hydrogel",
    "Microbes_Liquid" = "Microbes\nLiquid",
    "Neg_ctrl"        = "Neg. ctrl"
  )) +
  labs(
    title = "Wheat biomass partitioning by treatment",
    x     = NULL,
    y     = "Mean biomass (g)",
    fill  = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  )
biomasspartion_plot

ggsave("yield/plots/wheatbiomasspartioning.png", plot = biomasspartion_plot)

# ****************************************************
# Grain weight ANOVA ----
# ****************************************************

## with 8 tx data set ----
anova_grainweight8 <- aov(GrainWeight_g ~ Treatment, data = wheatyield_2025)
summary(anova_grainweight8)

tukey_grainweight8 <- TukeyHSD(anova_grainweight8)

# Shapiro Wilk test for normality
shapiro.test(anova_grainweight7$residuals)
# W = 0.97705, p-value = 0.4486 - data normally distributed

# Residuals distribution
hist(anova_grainweight7$residuals, 
     main = "Grain Weight Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_grainweight7$residuals, main = "Normal Q-Q Plot - Grain Weight Residuals")
qqline(anova_grainweight7$residuals, col = "lightblue", lwd = 2)

##12 tx data set ----
master <- read.csv("2025greenhousemaster_12tx.csv")

anova_grainweight <- aov(Grain_weight ~ Treatment_Category, data = master)
summary(anova_grainweight)

tukey_grainweight <- TukeyHSD(anova_grainweight)

# Shapiro Wilk test for normality
shapiro.test(anova_grainweight$residuals)
# W = 0.9729, p-value = 0.07332 - data normally distributed

# Residuals distribution
hist(anova_grainweight$residuals, 
     main = "Grain Weight Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_grainweight$residuals, main = "Normal Q-Q Plot - Grain Weight Residuals")
qqline(anova_grainweight$residuals, col = "lightblue", lwd = 2)

# ****************************************************
# Harvest Index ANOVA ----
# ****************************************************
## with 8 tx data set ----
anova_harvestindex8 <- aov(Harvest_Index ~ Treatment, data = wheatyield_2025)
summary(anova_harvestindex8)

tukey_harvestindex8 <- TukeyHSD(anova_harvestindex8)
tukey_harvestindex8

# Shapiro Wilk test for normality
shapiro.test(anova_harvestindex8$residuals)

# Residuals distribution
hist(anova_harvestindex8$residuals, 
     main = "Harvest Index Residuals", 
     xlab = "Residuals",
     col = "lightblue",
     breaks = 20)

# Quantile-Quantile plots
qqnorm(anova_harvestindex8$residuals, main = "Normal Q-Q Plot - Harvest Index Residuals")
qqline(anova_harvestindex8$residuals, col = "lightblue", lwd = 2)

# ****************************************************
# Combined bar graphs ----
# ****************************************************
AMF_yield <- read.csv("2025greenhousedata_focusedtx.csv", header = TRUE)

# summarize means and SE for all 4 variables
multi_summary <- AMF_yield %>%
  group_by(Treatment_Category) %>%
  summarise(
    mean_grain  = mean(GrainWeight_g,    na.rm = TRUE),
    se_grain    = sd(GrainWeight_g,      na.rm = TRUE) / sqrt(sum(!is.na(GrainWeight_g))),
    mean_hi     = mean(Harvest_Index,    na.rm = TRUE),
    se_hi       = sd(Harvest_Index,      na.rm = TRUE) / sqrt(sum(!is.na(Harvest_Index))),
    mean_rs     = mean(Root_Shoot,       na.rm = TRUE),
    se_rs       = sd(Root_Shoot,         na.rm = TRUE) / sqrt(sum(!is.na(Root_Shoot))),
    mean_amf    = mean(AMF_colonization, na.rm = TRUE),
    se_amf      = sd(AMF_colonization,   na.rm = TRUE) / sqrt(sum(!is.na(AMF_colonization))),
    .groups = "drop"
  )

# normalize means to max 
minmax <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

multi_norm <- multi_summary %>%
  mutate(
    norm_grain = minmax(mean_grain),
    norm_hi    = minmax(mean_hi),
    norm_rs    = minmax(mean_rs),
    norm_amf   = minmax(mean_amf)
  )

# normalize SE proportionally to the same scale as each variable
# SE is scaled by the same range used for each variable's min-max
multi_norm <- multi_norm %>%
  mutate(
    range_grain = max(mean_grain, na.rm=TRUE) - min(mean_grain, na.rm=TRUE),
    range_hi    = max(mean_hi,    na.rm=TRUE) - min(mean_hi,    na.rm=TRUE),
    range_rs    = max(mean_rs,    na.rm=TRUE) - min(mean_rs,    na.rm=TRUE),
    range_amf   = max(mean_amf,   na.rm=TRUE) - min(mean_amf,   na.rm=TRUE),
    se_norm_grain = se_grain / range_grain,
    se_norm_hi    = se_hi    / range_hi,
    se_norm_rs    = se_rs    / range_rs,
    se_norm_amf   = se_amf   / range_amf
  )

# pivot to long 
plot_data <- multi_norm %>%
  dplyr::select(Treatment, norm_grain, norm_hi, norm_rs, norm_amf,
                se_norm_grain, se_norm_hi, se_norm_rs, se_norm_amf) %>%
  tidyr::pivot_longer(
    cols      = c(norm_grain, norm_hi, norm_rs, norm_amf),
    names_to  = "Variable",
    values_to = "norm_mean"
  ) %>%
  mutate(
    norm_se = case_when(
      Variable == "norm_grain" ~ multi_norm$se_norm_grain[match(Treatment, multi_norm$Treatment)],
      Variable == "norm_hi"    ~ multi_norm$se_norm_hi[match(Treatment, multi_norm$Treatment)],
      Variable == "norm_rs"    ~ multi_norm$se_norm_rs[match(Treatment, multi_norm$Treatment)],
      Variable == "norm_amf"   ~ multi_norm$se_norm_amf[match(Treatment, multi_norm$Treatment)]
    ),
    Variable = dplyr::recode(Variable,
                             "norm_grain" = "Grain weight",
                             "norm_hi"    = "Harvest index",
                             "norm_rs"    = "Root:shoot ratio",
                             "norm_amf"   = "AMF colonization"
    ),
    Variable = factor(Variable, levels = c("Grain weight", "Harvest index",
                                           "Root:shoot ratio", "AMF colonization"))
  )

# Treatment order
treatment_order <- c("BioBead", "BioBead_SL", "Fertilizer",
                     "HardyGro", "Hydrogel", "Liquid_Microbes", "Negative")

plot_data <- plot_data %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order))

# Color palette 
variable_colors <- c(
  "Grain weight"     = "#E8834E",   
  "Harvest index"    = "#8DD3C7",   
  "Root:shoot ratio" = "#8B5E3C",   
  "AMF colonization" = "#BEBADA"    
)

# Plot
ggplot(plot_data, aes(x = Treatment, y = norm_mean, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = pmax(norm_mean - norm_se, 0),   
        ymax = norm_mean + norm_se),
    position = position_dodge(width = 0.75),
    width = 0.2, linewidth = 0.5, color = "grey30"
  ) +
  scale_fill_manual(values = variable_colors) +
  scale_x_discrete(labels = c(
    "BioBead"        = "BioBead",
    "BioBead_SL"     = "BioBead SL",
    "Fertilizer"     = "Fertilizer",
    "HardyGro"       = "HardyGro",
    "Hydrogel"       = "Hydrogel",
    "Liquid_Microbes"= "Microbes\nLiquid",
    "Negative"       = "Neg. ctrl"
  )) +
  scale_y_continuous(
    name   = "Min-Max normalized value (0–1)",
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Wheat response variables by treatment",
    x     = NULL,
    fill  = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(size = 9),
    plot.title      = element_text(hjust = 0.5, face = "bold")
  )

