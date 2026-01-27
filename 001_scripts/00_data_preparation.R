
# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: Part one of preparing data and general settings
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2025 06


# Global settings ---------------------------------------------------------
# Here global settings which are used throughout the script are defined

# Maximum lifespans (years). Taken from fishbase. 
age_max_AC <- 25
age_max_AS <- 54
age_max_EH <- 20
age_max_JM <- 5
age_max_ZF <- 5.5

# Settings ----------------------------------------------------------------
library(tidyverse)

# Colors -------------------------------------------------------------
# This section defines color palettes for the entire work
# Colors were generated using iwanthue online

# Colors for each species. 
color_species <- c("AC" = "#332288",
                   "AS" = "#DDCC77", 
                   "EH" = "#44AA99",
                   "JM" = "#88CCEE",
                   "ZF" = "#882255")
# For common name
color_species_common <- c("Atlantic cod" = "#332288",
                          "Australasian snapper" = "#DDCC77", 
                          "European hake" = "#44AA99",
                          "Japanese medaka" = "#88CCEE",
                          "Zebrafish" = "#882255")

# For scientific name
color_species_sci <- c("Gadus morhua" = "#332288",
                       "Chrysophrys auratus" = "#DDCC77", 
                       "Merluccius merluccius" = "#44AA99",
                       "Oryzias latipes" = "#88CCEE",
                       "Danio rerio" = "#882255")


# Colors to make direct comparisons between A and B 
color_compare <- c("#005AB5", "#DC3220")

# Color for sex
color_sex <- c("Females" = "#FFDA66", 
               "Males" = "#8ACE7E", 
               "Unknown" = "grey")

# Saving palettes 
save(color_compare, color_species, color_sex, file = "000_data/000_metadata/color_palettes.RData")

# Loading data ------------------------------------------------------------

# Loading prepared metadata for every species
# Atlantic cod (AC) - Gadus morhua
AC <- read_csv("000_data/000_metadata/AC_age.csv")

# Australasian snapper (AS) - Chrysophrys auratus
AS <- read_csv("000_data/000_metadata/AS_age.csv")

# European hake (EH) - Merluccius merluccius
EH <- read_csv("000_data/000_metadata/EH_age.csv")

# Japanese medaka (JM) - Oryzias latipes
JM <- read_csv("000_data/000_metadata/JM_age.csv")

# Zebrafish (ZF) - Danio rerio
ZF <- read_csv("000_data/000_metadata/ZF_age.csv")

# Data pre-processing ---------------------------------------------------------

# Creating one tibble with all species
meta_data <- rbind(AC, AS, EH, JM, ZF) %>% 
  # removing empty column
  select(-...1) %>%
  # Adding factors with levels
  mutate(name = factor(name, levels = sort(name)),
         sex = factor(str_to_title(sex), levels = c("Male", "Female", "Unknown")), 
         species_f = factor(species, levels = sort(unique(species))), .after = species) %>% 
  # Adding common name
  mutate(common_name = case_when(
    species == "AC" ~ "Atlantic cod",
    species == "AS" ~ "Australasian snapper", 
    species == "EH" ~ "European hake", 
    species == "JM" ~ "Japanese medaka",
    species == "ZF" ~ "Zebrafish")) %>% 
  # Getting rid of old columns
  select(-c(max_age, rel_age))


# Changing maximum and relative age
meta_data <- meta_data %>% 
  mutate(age_max = case_when(
    species == "AC" ~ age_max_AC,
    species == "AS" ~ age_max_AS, 
    species == "EH" ~ age_max_EH, 
    species == "JM" ~ age_max_JM,
    species == "ZF" ~ age_max_ZF)) %>% 
  relocate(age_max, .after = age) %>% 
  mutate(age_rel = age/age_max, .after = age_max) %>% 
  relocate(age_rel, .after = age_max)

# Saving metadata
save(meta_data, file = "000_data/000_metadata/metadata.RData")
