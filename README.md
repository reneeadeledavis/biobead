# Introduction

BioBead is a biodegradable hydrogel pellet containing fungal and
bacterial consortia for soil health and regenerative agriculture.
Microbes in this consortia can fix nitrogen (N), solubilize phosphorus
(P), and drive carbon (C) cycling. Arbuscular mycorrhizal fungi (AMF)
are obligate mutualists that make direct contact with plant roots and
hand off N and P to plants, obtaining C in return. They also provide
structures for microbial residence and travel. Greenhouse studies with
BioBead have demonstrated improved growth for tomatoes and wheat plants,
increases in soil organic matter, and promotion of closed nutrient
cycling.

In 2025, a greenhouse study led by Renee Davis and Korena Mafune, PhD
evaluated the impacts of 14 different microbe inoculations and
treatments on soil microbiology, wheat growth, root colonization, and
greenhouse gas emissions. This specific project will focus on soil
geochemistry, nutrient dynamics, and soil moisture data, with other data
being incorporated in subsequent phases. The responses measured at
different intervals will have to be analyzed separately, or combined in
another manner. These results will provide a foundation for expansion
into larger field trials and BioBead formulation improvements in 2026
and beyond, and hopefully become publications and dissertation chapters.

# Methods

## Greenhouse Trial and Sample Collection

Wheat seedlings (Red Fife variety, Territorial Seeds) were planted in
1-gallon pots containing sunshine mix and sand (1:2 v/v) from
September-December 2025 at the UW Biology Greenhouse. The experiment
utilized a randomized design with 14 treatments and 7 biological
replicates per treatment (98 pots total). Microbial inoculants were
applied at planting according to treatment specifications.

Plants were grown for 14 weeks under greenhouse conditions. Plants were
watered weekly with 200 mL water, and soil moisture was recorded 24
hours after each watering event. Leachate samples were collected at
weeks 0, 2, 4, 6, 8 and analyzed using The Gallery (ThermoFisher) for
ammonium, nitrate, and phosphate concentrations.

Plant Root Simulator (PRS) Probes (Western Ag Innovations) were buried
at week 6 and retrieved after 14 days. Probes were sent to Western Ag
for analysis of soil cation and anion exchange capacity via HPLC. At
experiment termination (week 14), above-ground biomass was harvested and
weighed, and rhizosphere soil and root samples were collected and stored
at -80°C for future DNA extraction and metagenomic sequencing. (These
data are mentioned for reference and will not be analyzed in this
project.)

# Data

## Files

**2025biobeadgreenhouse.Rmd** - R Markdown file containing this
document (project overview). 

**2025greenhousemaster.csv** - master data file for multivariate statistical analysis.

**Leachates_2025.csv** - Ammonium, nitrate, and phosphate leachate
concentrations recorded throughout the study. (This is a time series, so
temporal autocorrelation may be important.) There is an associated .Rmd file with an exploratory one-way ANOVA. 

**PRS_GH_2025.csv** - Soil anion and cation exchange data from PRS
probes. There is an associated .Rmd file with an exploratory one-way ANOVA. 

**WheatMoistureData_2025.csv** - Soil moisture data recorded throughout
the study.(This is a time series, so temporal autocorrelation may be
important.) There is an associated .Rmd file with an exploratory one-way ANOVA. 

**wheat yield.csv** - Wheat yield data from end of experiment. 

**N2O_fluxdata.csv** - Nitrous oxide flux reported as N"[2]*"O-N/sqmt/hr. Calculated from peak area. 

## Response variables

| Code             | Full Name             | Units              | Range      |
|------------------|-----------------------|--------------------|------------|
| NH4-N_ugml       | Ammonium              | ug/mL              | 0.01-5.24  |
| NO3-N_ugml       | Nitrate               | ug/mL              | 0-97.55    |
| PO4-P_ugml       | Phosphate             | ug/mL              | 0.05-31.97 |
| \% soil moisture | Percent Soil Moisture | n/a                | 0-100      |
| NO3-N            | Nitrate               | ug/10cm^2^/14 days | 0.8-8.3    |
| NH4-N            | Ammonium              | ug/10cm^2^/14 days | 0-2.9      |
| P                | Phosphorus            | ug/10cm^2^/14 days | 6.8-23.7   |
| Ca               | Calcium               | ug/10cm^2^/14 days | 751-2220   |
| Mg               | Magnesium             | ug/10cm^2^/14 days | 158-477    |
| K                | Potassium             | ug/10cm^2^/14 days | 146-305    |
| Fe               | Iron                  | ug/10cm^2^/14 days | 2.5-21.5   |
| Mn               | Manganese             | ug/10cm^2^/14 days | 0.3-1.5    |
| Cu               | Copper                | ug/10cm^2^/14 days | 0-1        |
| Zn               | Zinc                  | ug/10cm^2^/14 days | 1.3-3.8    |
| B                | Boron                 | ug/10cm^2^/14 days | 0-0.6      |
| S                | Sulphur               | ug/10cm^2^/14 days | 105-490    |
| Pb               | Lead                  | ug/10cm^2^/14 days | 0-0.9      |
| Al               | Aluminum              | ug/10cm^2^/14 days | 3.4-9.5    |
| Cd               | Cadmium               | ug/10cm^2^/14 days | 0-0.2      |
| Na               | Sodium                | ug/10cm^2^/14 days | 943-2734   |
| GrainWeight_g | Grain weight         | g                | 0.6-2.7  |
| Harvest_Index | Harvest index        | n/a              | 0-0.6    |
| N2O-N_flux    | Nitrous oxide N flux | ug N2O-N/sqmt/hr | -38-3557 |

## Explanatory variables

| Code | Full Name | Units | Range/Values |
|----|----|----|----|
| Fungi | Arbuscular mycorrhizal fungi added | Binary | Present/absent |
| Bacteria | Bacterial consortia added | Binary | Present/absent |
| Hydrogel | BioBead hydrogel carrier | Binary | Present/absent |
| Dose | Inoculum dose | Categorical | Full (100%), Half (50%) |
| GR24 | Synthetic strigolactone addition (5 nM) | Binary | Present/absent |
| Week | Sampling week | Integer | 1-8 |
| Treatment | Treatment code | Categorical | T1-T14 (see below) |

## Treatment codes

T1 - BioBead – bacteria

T2 - BioBead – AMF

T3 - BioBead - bacteria, AMF

T4 - SL-primed BioBead - bacteria, AMF

T5 - SL-primed BioBead - bacteria

T6 - SL-primed BioBead - AMF

T7 - BioBead w/ reapplication at week 8

T8 - ½ fertilizer ½ BioBead

T9 - Liquid consortia

T10 - Commercial biofertilizer (Hardygro)

T11 - Empty bead

T12 - SL-primed empty bead

T13 - Synthetic fertilizer (positive control – 10% dose, 0.6 g/L,
10mL/pot)

T14 - Untreated

Specific replicates of samples are indicated after a dash in the sample
code. For example, certain replicates in the PRS data are designated as
T1-1, T1-2, T1-3, etc.
