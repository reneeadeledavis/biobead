# ****************************************************
# Load packages ----
# ****************************************************
library(vegan)
library(dplyr)

# ****************************************************
# Read in data ----
# ****************************************************
master <- read.csv("2025greenhousemaster_12tx.csv")

# ****************************************************
# 2 factor PERMANOVA ----
# ****************************************************

## SL priming effects on AMF and plant responses ----
# Filter to bead treatments only to balance factorials
bead_treatments <- c("BioBead_bac", "BioBead_AMF", "BioBead_Mixed",
                     "Empty_Bead", "SL_BioBead_Bac", "SL_BioBead_AMF",
                     "SL_BioBead_Mixed", "SL_Empty_Bead")

master_beads <- master %>%
  filter(Treatment_Category %in% bead_treatments) %>%
  mutate(
    SL  = ifelse(grepl("^SL", Treatment_Category), "Y", "N"),
    AMF = ifelse(grepl("AMF|Mixed", Treatment_Category), "Y", "N")
  )
head(master_beads) #verify coding is correct

table(master_beads$SL, master_beads$AMF) # Verify the factorial structure

properties <- master_beads[, c("AMF_colonization", "Arb_colonization",
                               "Grain_weight", "Root_biomass", "Total_biomass")]
properties_rel <- decostand(properties, method = "max")
data_dist <- vegdist(properties_rel, method = "euclidean")

permanova_2factor <- adonis2(data_dist ~ SL * AMF, 
                             data = master_beads,
                             permutations = 999,
                             by = "terms")
print(permanova_2factor)

## SL priming effects on AMF establishment ----
amf_properties <- master_beads[, c("Hyph_colonization", "Arb_colonization")]
amf_rel <- decostand(amf_properties, method = "max")
amf_dist <- vegdist(amf_rel, method = "euclidean")

permanova_amf <- adonis2(amf_dist ~ SL * AMF,
                         data = master_beads,
                         permutations = 999,
                         by = "terms")
print(permanova_amf)

## SL priming effects on plant responses ----
plant_properties <- master_beads[, c("Grain_weight", "Root_biomass", "Total_biomass")]
plant_rel <- decostand(plant_properties, method = "max")
plant_dist <- vegdist(plant_rel, method = "euclidean")

permanova_plant <- adonis2(plant_dist ~ SL * AMF,
                           data = master_beads,
                           permutations = 999,
                           by = "terms")
print(permanova_plant)

## SL priming effects on whole system ----
whole_properties <- master_beads[, c("Hyph_colonization", "Arb_colonization",
                                     "Grain_weight", "Root_biomass", "Total_biomass",
                                     "NH4.N_ugml", "NO3.N_ugml", "PO4.P_ugml")]
whole_rel <- decostand(whole_properties, method = "max")
whole_dist <- vegdist(whole_rel, method = "euclidean")

permanova_whole <- adonis2(whole_dist ~ SL * AMF,
                           data = master_beads,
                           permutations = 999,
                           by = "terms")
print(permanova_whole)

# ****************************************************
# 2 factor ANOVA ----
# ****************************************************

# SL and AMF colonization
anova_sl_amf <- aov(AMF_colonization ~ SL * AMF, data = master_beads)
summary(anova_sl_amf)

# SL and arbuscule formation 
anova_sl_arb <- aov(Arb_colonization ~ SL * AMF, data = master_beads)
summary(anova_sl_arb)

# SL and grain weight
anova_sl_grain <- aov(Grain_weight ~ SL * AMF, data = master_beads)
summary(anova_sl_grain)

# ****************************************************
# Single factor ANOVA ----
# ****************************************************

# note: the objects are named identical to the section above, so running these commands will override them
# SL and AMF colonization 
anova_sl_amf <- aov(AMF_colonization ~ SL, data = master_beads)
summary(anova_sl_amf)

# SL and intraradical hyphae
anova_sl_hyph <- aov(Hyph_colonization ~ SL, data = master_beads)
summary(anova_sl_hyph)

# SL and arbuscule formation 
anova_sl_arb <- aov(Arb_colonization ~ SL, data = master_beads)
summary(anova_sl_arb)

# SL and grain weight
anova_sl_grain <- aov(Grain_weight ~ SL, data = master_beads)
summary(anova_sl_grain)
