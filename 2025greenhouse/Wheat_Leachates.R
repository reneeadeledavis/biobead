##Load packages
library(ggstatsplot)
library(ggplot2)
library(reshape)
library(dplyr)
library(plotrix)
library(see)
library(ggpubr)

##Read in data
Leachates_2025 <- read.csv("Leachates_2025.csv", header = TRUE)

#Summary stats
Summary_all_weeks <- Leachates_2025 %>% 
  group_by(Date, Treatment_Category) %>%
  summarize(
    NH4_mean = mean(NH4.N_ugml, na.rm = TRUE),
    NH4_sd = sd(NH4.N_ugml, na.rm = TRUE),
    NH4_se = sd(NH4.N_ugml, na.rm = TRUE) / sqrt(sum(!is.na(NH4.N_ugml))),
    NO3_mean = mean(NO3.N_ugml, na.rm = TRUE),
    NO3_sd = sd(NO3.N_ugml, na.rm = TRUE),
    NO3_se = sd(NO3.N_ugml, na.rm = TRUE) / sqrt(sum(!is.na(NO3.N_ugml))),
    PO4_mean = mean(PO4.P_ugml, na.rm = TRUE),
    PO4_sd = sd(PO4.P_ugml, na.rm = TRUE),
    PO4_se = sd(PO4.P_ugml, na.rm = TRUE) / sqrt(sum(!is.na(PO4.P_ugml))),
    n = n(),
    .groups = 'drop'
  )

#Palette for plots
JACE_COLOR_2 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                  "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", 
                  "#bcbd22", "#17becf", )

#Subset data by week
Week2_subset <- subset(Leachates_2025, Date=="Week_2")
Week4_subset <- subset(Leachates_2025, Date=="Week_4")
Week6_subset <- subset(Leachates_2025, Date=="Week_6")
Week8_subset <- subset(Leachates_2025, Date=="Week_8")

# ****************************************************
# Week 2 ----
# ****************************************************

###Data transformation/melting
melted_leachates2 <- melt(Week2_subset , 
                         id.vars = c("Date", "Sample_ID", "Treatment_Category"), 
                         measure.vars = c("NH4.N_ugml", "NO3.N_ugml", "PO4.P_ugml"),
                         variable.name = "Extractant", 
                         value.name = "Concentration")
write.csv(melted_leachates2, "Week2_melt_2025.csv")
wk2_melt <- read.csv("Week2_melt_2025.csv", header = TRUE) #seems redundant with line 17?

Leachates_2025$NO3.N_ug <- as.numeric(as.character(Leachates_2025$NO3.N_ug)) #not needed given how R is reading NO3 data as numeric

##Generating stats
Summary_leachateswk2 <- Week2_subset %>%
  group_by(Treatment_Category) %>%
  summarize(NH4_mean=mean(NH4.N_ugml, na.rm=TRUE), 
            NO3_mean=mean(NO3.N_ugml, na.rm=TRUE), 
            PO4_mean=mean(PO4.P_ugml, na.rm=TRUE), 
            NH4_se=std.error(NH4.N_ugml, na.rm=TRUE),
            NO3_se=std.error(NO3.N_ugml, na.rm=TRUE),
            PO4_se=std.error(PO4.P_ugml, na.rm=TRUE),
            NH4_max=max(NH4.N_ugml, na.rm=TRUE), 
            NO3_max=max(NO3.N_ugml, na.rm=TRUE), 
            PO4_max=max(PO4.P_ugml, na.rm=TRUE), 
            NH4_min=min(NH4.N_ugml, na.rm=TRUE), 
            NO3_min=min(NO3.N_ugml, na.rm=TRUE), 
            PO4_min=min(PO4.P_ugml, na.rm=TRUE))

###Violin plot
ggplot(wk2_melt, aes(x = Treatment_Category, y = value, Group = variable))+
  geom_violinhalf(position = position_nudge(x = .2), 
                  alpha = .4) +
  geom_point(aes(color = variable), 
             position = position_jitter(width = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  scale_fill_manual(values = JACE_COLOR_2) +
  scale_color_manual(values = JACE_COLOR_2) +
  labs(title="Week 2 Leachate", x="Treatment", y = "Concentration (mg/ml)") +
  theme_bw() +
  theme(legend.position="none") +
  geom_jitter(shape=20, position=position_jitter(0.1)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.75, hjust=0.8))+
  theme(plot.title = element_text(hjust = 0.5, size = 22)) +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

###Faceted violin plot
ggplot(wk2_melt, aes(x = Treatment_Category, y = value, fill = Treatment_Category)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_jitter(width = 0.1, 
              size = 1.5, 
              alpha = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Horizontal layout
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Week 2 Nutrient Leachates", 
       x = "Treatment", 
       y = "Concentration (µg/mL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

###Correlation scatterplots
ggscatter(Week2_subset, x = "PO4.P_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PO4", ylab = "NO3")

ggscatter(Week2_subset, x = "NH4.N_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NH4", ylab = "NO3")

# ****************************************************
# Week 4 ----
# ****************************************************

###Data transformation/melting
melted_leachates4 <- melt(Week4_subset , 
                          id.vars = c("Date", "Sample_ID", "Treatment_Category"), 
                          measure.vars = c("NH4.N_ugml", "NO3.N_ugml", "PO4.P_ugml"),
                          variable.name = "Extractant", 
                          value.name = "Concentration")
write.csv(melted_leachates4, "Week4_melt_2025.csv")
wk4_melt <- read.csv("Week4_melt_2025.csv", header = TRUE) #seems redundant with line 17?

Leachates_2025$NO3.N_ug <- as.numeric(as.character(Leachates_2025$NO3.N_ug)) #not needed given how R is reading NO3 data as numeric

##Generating stats
Summary_leachateswk4 <- Week4_subset %>%
  group_by(Treatment_Category) %>%
  summarize(NH4_mean=mean(NH4.N_ugml, na.rm=TRUE), 
            NO3_mean=mean(NO3.N_ugml, na.rm=TRUE), 
            PO4_mean=mean(PO4.P_ugml, na.rm=TRUE), 
            NH4_se=std.error(NH4.N_ugml, na.rm=TRUE),
            NO3_se=std.error(NO3.N_ugml, na.rm=TRUE),
            PO4_se=std.error(PO4.P_ugml, na.rm=TRUE),
            NH4_max=max(NH4.N_ugml, na.rm=TRUE), 
            NO3_max=max(NO3.N_ugml, na.rm=TRUE), 
            PO4_max=max(PO4.P_ugml, na.rm=TRUE), 
            NH4_min=min(NH4.N_ugml, na.rm=TRUE), 
            NO3_min=min(NO3.N_ugml, na.rm=TRUE), 
            PO4_min=min(PO4.P_ugml, na.rm=TRUE))

###Violin plot
ggplot(wk4_melt, aes(x = Treatment_Category, y = value, Group = variable))+
  geom_violinhalf(position = position_nudge(x = .2), 
                  alpha = .4) +
  geom_point(aes(color = variable), 
             position = position_jitter(width = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  scale_fill_manual(values = JACE_COLOR_2) +
  scale_color_manual(values = JACE_COLOR_2) +
  labs(title="Week 4 Leachate", x="Treatment", y = "Concentration (mg/ml)") +
  theme_bw() +
  theme(legend.position="none") +
  geom_jitter(shape=20, position=position_jitter(0.1)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.75, hjust=0.8))+
  theme(plot.title = element_text(hjust = 0.5, size = 22)) +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

###Faceted violin plot
ggplot(wk4_melt, aes(x = Treatment_Category, y = value, fill = Treatment_Category)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_jitter(width = 0.1, 
              size = 1.5, 
              alpha = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Horizontal layout
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Week 4 Nutrient Leachates", 
       x = "Treatment", 
       y = "Concentration (µg/mL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

###Correlation scatterplots
ggscatter(Week4_subset, x = "PO4.P_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PO4", ylab = "NO3")

ggscatter(Week4_subset, x = "NH4.N_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NH4", ylab = "NO3")

# ****************************************************
# Week 6 ----
# ****************************************************

###Data transformation/melting
melted_leachates6 <- melt(Week6_subset , 
                          id.vars = c("Date", "Sample_ID", "Treatment_Category"), 
                          measure.vars = c("NH4.N_ugml", "NO3.N_ugml", "PO4.P_ugml"),
                          variable.name = "Extractant", 
                          value.name = "Concentration")
write.csv(melted_leachates6, "Week6_melt_2025.csv")
wk6_melt <- read.csv("Week6_melt_2025.csv", header = TRUE) #seems redundant with line 17?

Leachates_2025$NO3.N_ug <- as.numeric(as.character(Leachates_2025$NO3.N_ug)) #not needed given how R is reading NO3 data as numeric

##Generating stats
Summary_leachateswk6 <- Week6_subset %>%
  group_by(Treatment_Category) %>%
  summarize(NH4_mean=mean(NH4.N_ugml, na.rm=TRUE), 
            NO3_mean=mean(NO3.N_ugml, na.rm=TRUE), 
            PO4_mean=mean(PO4.P_ugml, na.rm=TRUE), 
            NH4_se=std.error(NH4.N_ugml, na.rm=TRUE),
            NO3_se=std.error(NO3.N_ugml, na.rm=TRUE),
            PO4_se=std.error(PO4.P_ugml, na.rm=TRUE),
            NH4_max=max(NH4.N_ugml, na.rm=TRUE), 
            NO3_max=max(NO3.N_ugml, na.rm=TRUE), 
            PO4_max=max(PO4.P_ugml, na.rm=TRUE), 
            NH4_min=min(NH4.N_ugml, na.rm=TRUE), 
            NO3_min=min(NO3.N_ugml, na.rm=TRUE), 
            PO4_min=min(PO4.P_ugml, na.rm=TRUE))

###Violin plot
ggplot(wk6_melt, aes(x = Treatment_Category, y = value, Group = variable))+
  geom_violinhalf(position = position_nudge(x = .2), 
                  alpha = .4) +
  geom_point(aes(color = variable), 
             position = position_jitter(width = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  scale_fill_manual(values = JACE_COLOR_2) +
  scale_color_manual(values = JACE_COLOR_2) +
  labs(title="Week 6 Leachate", x="Treatment", y = "Concentration (mg/ml)") +
  theme_bw() +
  theme(legend.position="none") +
  geom_jitter(shape=20, position=position_jitter(0.1)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.75, hjust=0.8))+
  theme(plot.title = element_text(hjust = 0.5, size = 22)) +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

###Faceted violin plot
ggplot(wk6_melt, aes(x = Treatment_Category, y = value, fill = Treatment_Category)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_jitter(width = 0.1, 
              size = 1.5, 
              alpha = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Horizontal layout
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Week 6 Nutrient Leachates", 
       x = "Treatment", 
       y = "Concentration (µg/mL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

###Correlation scatterplots
ggscatter(Week6_subset, x = "PO4.P_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PO4", ylab = "NO3")

ggscatter(Week6_subset, x = "NH4.N_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NH4", ylab = "NO3")

# ****************************************************
# Week 8 ----
# ****************************************************

###Data transformation/melting
melted_leachates8 <- melt(Week8_subset , 
                          id.vars = c("Date", "Sample_ID", "Treatment_Category"), 
                          measure.vars = c("NH4.N_ugml", "NO3.N_ugml", "PO4.P_ugml"),
                          variable.name = "Extractant", 
                          value.name = "Concentration")
write.csv(melted_leachates8, "Week8_melt_2025.csv")
wk8_melt <- read.csv("Week8_melt_2025.csv", header = TRUE) #seems redundant with line 17?

Leachates_2025$NO3.N_ug <- as.numeric(as.character(Leachates_2025$NO3.N_ug)) #not needed given how R is reading NO3 data as numeric

##Generating stats
Summary_leachateswk8 <- Week8_subset %>%
  group_by(Treatment_Category) %>%
  summarize(NH4_mean=mean(NH4.N_ugml, na.rm=TRUE), 
            NO3_mean=mean(NO3.N_ugml, na.rm=TRUE), 
            PO4_mean=mean(PO4.P_ugml, na.rm=TRUE), 
            NH4_se=std.error(NH4.N_ugml, na.rm=TRUE),
            NO3_se=std.error(NO3.N_ugml, na.rm=TRUE),
            PO4_se=std.error(PO4.P_ugml, na.rm=TRUE),
            NH4_max=max(NH4.N_ugml, na.rm=TRUE), 
            NO3_max=max(NO3.N_ugml, na.rm=TRUE), 
            PO4_max=max(PO4.P_ugml, na.rm=TRUE), 
            NH4_min=min(NH4.N_ugml, na.rm=TRUE), 
            NO3_min=min(NO3.N_ugml, na.rm=TRUE), 
            PO4_min=min(PO4.P_ugml, na.rm=TRUE))

###Violin plot
ggplot(wk8_melt, aes(x = Treatment_Category, y = value, Group = variable))+
  geom_violinhalf(position = position_nudge(x = .2), 
                  alpha = .4) +
  geom_point(aes(color = variable), 
             position = position_jitter(width = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5) +
  scale_fill_manual(values = JACE_COLOR_2) +
  scale_color_manual(values = JACE_COLOR_2) +
  labs(title="Week 8 Leachate", x="Treatment", y = "Concentration (mg/ml)") +
  theme_bw() +
  theme(legend.position="none") +
  geom_jitter(shape=20, position=position_jitter(0.1)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.75, hjust=0.8))+
  theme(plot.title = element_text(hjust = 0.5, size = 22)) +
  theme(axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) 

###Faceted violin plot
ggplot(wk8_melt, aes(x = Treatment_Category, y = value, fill = Treatment_Category)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.2, 
               outlier.shape = NA,
               alpha = 0.5) +
  geom_jitter(width = 0.1, 
              size = 1.5, 
              alpha = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Horizontal layout
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Week 8 Nutrient Leachates", 
       x = "Treatment", 
       y = "Concentration (µg/mL)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

###Correlation scatterplots
ggscatter(Week8_subset, x = "PO4.P_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "PO4", ylab = "NO3")

ggscatter(Week8_subset, x = "NH4.N_ugml", y = "NO3.N_ugml", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NH4", ylab = "NO3")
