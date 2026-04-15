library(dplyr)
library(tidyr)
library(ggplot2)

master <- read.csv("2025greenhousemaster_8tx.csv")

treatment_colors <- c(
  "BioBead_AMF"    = "#8DD3C7",
  "BioBead_Bac"    = "#BC80BD",
  "BioBead_Mixed"  = "#BEBADA",
  "Fertilizer"     = "#FB8072",
  "Liquid_Microbes"= "#FFFFB3",
  "HardyGro"       = "#FDB462",
  "Empty_Bead"     = "#80B1D3",
  "Negative"       = "#D9D9D9"
)

leachate_summary <- master %>%
  group_by(Treatment_Category) %>%
  summarise(
    mean_NH4 = mean(`NH4.N_ugml`, na.rm = TRUE),
    mean_NO3 = mean(`NO3.N_ugml`, na.rm = TRUE),
    mean_PO4 = mean(`PO4.P_ugml`, na.rm = TRUE),
    se_NH4   = sd(`NH4.N_ugml`, na.rm = TRUE) / sqrt(sum(!is.na(`NH4.N_ugml`))),
    se_NO3   = sd(`NO3.N_ugml`, na.rm = TRUE) / sqrt(sum(!is.na(`NO3.N_ugml`))),
    se_PO4   = sd(`PO4.P_ugml`, na.rm = TRUE) / sqrt(sum(!is.na(`PO4.P_ugml`))),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = c(mean_NH4, mean_NO3, mean_PO4),
    names_to  = "Nutrient",
    values_to = "Mean_ugml"
  ) %>%
  mutate(
    SE = case_when(
      Nutrient == "mean_NH4" ~ se_NH4,
      Nutrient == "mean_NO3" ~ se_NO3,
      Nutrient == "mean_PO4" ~ se_PO4
    ),
    Nutrient = recode(Nutrient,
                      "mean_NH4" = "NH4-N",
                      "mean_NO3" = "NO3-N",
                      "mean_PO4" = "PO4-P"
    )
  )

ggplot(leachate_summary, aes(x = Treatment_Category, y = Mean_ugml, fill = Treatment_Category)) +
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
