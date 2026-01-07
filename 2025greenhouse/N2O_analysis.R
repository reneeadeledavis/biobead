# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read in data
data <- read_excel("GHG_wheat_2025_dataarrangement.xlsx")

# Display data structure
cat("=== Data Structure ===\n")
str(data)
cat("\n")

# Calculate summary statistics by treatment and sample period
summary_data <- data %>%
  group_by(Sample_Period, Treatment_Level) %>%
  summarise(
    mean_area = mean(Curve_Area, na.rm = TRUE),
    sd_area = sd(Curve_Area, na.rm = TRUE),
    se_area = sd(Curve_Area, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  )

# Display summary statistics
cat("=== Summary Statistics ===\n")
print(summary_data, n = Inf)
cat("\n")

# Create output directory 
if (!dir.exists("output")) {
  dir.create("output")
}

# 9/25 bar graph
cat("Creating September 25 plot...\n")

sept_data <- summary_data %>% 
  filter(Sample_Period == "Sept_25") %>%
  arrange(desc(mean_area))

p1 <- ggplot(sept_data, aes(x = reorder(Treatment_Level, -mean_area), y = mean_area)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_area - se_area, ymax = mean_area + se_area),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("N"[2]*"O Peak Area by Treatment - September 25, 2025")),
    x = "Treatment",
    y = expression(bold("N"[2]*"O Peak Area (Mean ± SE)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("output/N2O_Sept25_barplot.png", plot = p1, 
       width = 12, height = 7, dpi = 300, bg = "white")
cat("✓ September 25 plot saved to: output/N2O_Sept25_barplot.png\n\n")

# 10/11 bar graph
cat("Creating October 11 plot...\n")

oct_data <- summary_data %>% 
  filter(Sample_Period == "Oct_11") %>%
  arrange(desc(mean_area))

p2 <- ggplot(oct_data, aes(x = reorder(Treatment_Level, -mean_area), y = mean_area)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_area - se_area, ymax = mean_area + se_area),
                width = 0.3, linewidth = 0.5) +
  labs(
    title = expression(bold("N"[2]*"O Peak Area by Treatment - October 11, 2025")),
    x = "Treatment",
    y = expression(bold("N"[2]*"O Peak Area (Mean ± SE)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("output/N2O_Oct11_barplot.png", plot = p2, 
       width = 12, height = 7, dpi = 300, bg = "white")
cat("✓ October 11 plot saved to: output/N2O_Oct11_barplot.png\n\n")

# grouped bars
cat("Creating combined comparison plot...\n")

# Order treatments by average across both dates
treatment_order <- summary_data %>%
  group_by(Treatment_Level) %>%
  summarise(avg_mean = mean(mean_area)) %>%
  arrange(desc(avg_mean)) %>%
  pull(Treatment_Level)

summary_data_ordered <- summary_data %>%
  mutate(Treatment_Level = factor(Treatment_Level, levels = treatment_order))

p3 <- ggplot(summary_data_ordered, 
             aes(x = Treatment_Level, y = mean_area, fill = Sample_Period)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_area - se_area, ymax = mean_area + se_area),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.5) +
  scale_fill_manual(
    values = c("Sept_25" = "steelblue", "Oct_11" = "coral"),
    labels = c("Sept_25" = "September 25", "Oct_11" = "October 11"),
    name = "Collection Date"
  ) +
  labs(
    title = expression(bold("N"[2]*"O Peak Area Comparison - Both Collection Dates")),
    x = "Treatment",
    y = expression(bold("N"[2]*"O Peak Area (Mean ± SE)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 15, margin = margin(b = 15)),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("output/N2O_Combined_barplot.png", plot = p3, 
       width = 14, height = 8, dpi = 300, bg = "white")
cat("✓ Combined comparison plot saved to: output/N2O_Combined_barplot.png\n\n")

# facted plot
cat("Creating faceted plot...\n")

# Reorder within each facet
summary_data_facet <- summary_data %>%
  group_by(Sample_Period) %>%
  arrange(Sample_Period, desc(mean_area)) %>%
  mutate(order = row_number()) %>%
  ungroup()

p4 <- ggplot(summary_data_facet, 
             aes(x = reorder(Treatment_Level, -mean_area), y = mean_area)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7, 
           color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_area - se_area, ymax = mean_area + se_area),
                width = 0.3, linewidth = 0.5) +
  facet_wrap(~ Sample_Period, ncol = 1, scales = "free_x",
             labeller = labeller(Sample_Period = c(
               "Sept_25" = "September 25, 2025",
               "Oct_11" = "October 11, 2025"
             ))) +
  labs(
    title = expression(bold("N"[2]*"O Peak Area by Treatment - Both Collection Dates")),
    x = "Treatment",
    y = expression(bold("N"[2]*"O Peak Area (Mean ± SE)"))
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text= element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("output/N2O_Faceted_barplot.png", plot = p4, 
       width = 12, height = 10, dpi = 300, bg = "white")
cat("✓ Faceted plot saved to: output/N2O_Faceted_barplot.png\n\n")
