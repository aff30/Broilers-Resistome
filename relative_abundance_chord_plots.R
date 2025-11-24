library(phyloseq)
library(tidyverse)
library(circlize)
library(tibble)
library(forcats)
library(stringr)

ps <- readRDS("ps/decontam-ps.rds")
ps_amr <- ps

# Check sample data
sample_data(ps_amr)[1:3, ]
tax_table(ps_amr)[1:3, ]

# All AMR, no filtering yet
amr_long <- phyloseq::psmelt(ps_amr) %>%
  filter(Abundance > 0) %>%
  filter(!is.na(Treatment), !is.na(Class), !is.na(Age_Days))   # if your column is Age_Days, change here

amr_long <- amr_long %>%
  mutate(
    # 1) Capitalize Betalactams
    Class = ifelse(
      tolower(Class) == "betalactams",
      "Betalactams",
      Class
    ))
    

    
    amr_long <- amr_long %>%
      mutate(
        Class = str_wrap(Class, width = 20)  # automatically adds new lines
      )
    

BROADGROUP_COL <- "Broadclass"
DRUG_LABEL     <- "Drugs"

# ---- 1) Drug-only subset, with Age + Treatment -------------------------

amr_long_drug <- amr_long %>%
  filter(
    !is.na(.data[[BROADGROUP_COL]]),
    .data[[BROADGROUP_COL]] == DRUG_LABEL
  ) %>%
  # recode BMD -> Antibiotic
  mutate(
    Treatment = fct_recode(Treatment, "Antibiotic" = "BMD"),
    Treatment = as.character(Treatment),
    Days      = as.character(Age_Days)   # or Age_Days if that's your column
  )

# ---- 2) Summarise by Treatment × Day × Class ---------------------------

amr_sum_drug_age <- amr_long_drug %>%
  group_by(Treatment, Days, Class) %>%
  summarise(
    total_abund = sum(Abundance),
    .groups     = "drop"
  ) %>%
  mutate(
    weight  = log10(total_abund + 1),
    Trt_Day = paste(Treatment, paste0("D", Days), sep = "_")
    # e.g., "Basal Diet_D38"
  )

# ---- 3) Wide matrix: rows = Treatment_Day, cols = Drug AMR classes -----

amr_mat_drug_age <- amr_sum_drug_age %>%
  select(Trt_Day, Class, weight) %>%
  pivot_wider(
    names_from  = Class,
    values_from = weight,
    values_fill = 0
  ) %>%
  column_to_rownames("Trt_Day") %>%
  as.matrix()

dim(amr_mat_drug_age)
amr_mat_drug_age[, 1:min(5, ncol(amr_mat_drug_age))]

# ---- 4) Colors for Treatment×Day sectors and AMR classes ---------------

treat_colors <- c(
  "Basal Diet"     = "#eae69e",
  "Probiotic"      = "#bfdb81",
  "Essential oils" = "#83a561",
  "Antibiotic"     = "#48723e"
)

# Get all unique Trt_Day labels and map to base treatment colors
trt_day_labels <- rownames(amr_mat_drug_age)
trt_day_treat  <- sub("_D.*$", "", trt_day_labels)  # everything before "_D"

trt_day_colors <- setNames(
  treat_colors[trt_day_treat],
  trt_day_labels
)

# Drug class colors
amr_classes_drug <- colnames(amr_mat_drug_age)

amr_class_colors_drug <- setNames(
  grDevices::colorRampPalette(c("#6c757d", "#adb5bd", "#ced4da"))(length(amr_classes_drug)),
  amr_classes_drug
)

grid_colors_drug_age <- c(trt_day_colors, amr_class_colors_drug)

# ---- 5) Long links table for chordDiagram -------------------------------

links_drug_age <- amr_mat_drug_age %>%
  as.data.frame() %>%
  rownames_to_column("Trt_Day") %>%
  pivot_longer(
    cols      = -Trt_Day,
    names_to  = "Class",
    values_to = "weight"
  ) %>%
  filter(weight > 0)

head(links_drug_age)

# ---- 6) Chord diagram: Treatment×Age ↔ AMR drug classes ----------------

circos.clear()

chordDiagram(
  x               = links_drug_age,
  grid.col        = grid_colors_drug_age,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    circos.text(
      x      = CELL_META$xcenter,
      y      = CELL_META$ylim[1] - mm_y(4),
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj    = c(0, 0.5),
      cex    = 0.55  # smaller, more sectors
    )
  },
  bg.border = NA
)

title("Chord diagram: Treatment × Age \u2194 AMR drug classes (log10 abundance)")


##_______________________________________________________________________________##

## selecting the biocides and metals 

library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(circlize)

## ---- 1. Subset to Biocides + Metals -----------------------------------

BIO_MET_LABELS <- c("Biocides", "Metals")  # adjust if names differ in your taxonomy

amr_long_biomet <- amr_long %>%
  filter(!is.na(.data[[BROADGROUP_COL]]),
         .data[[BROADGROUP_COL]] %in% BIO_MET_LABELS,
         !is.na(Class))

## ---- 2. Summarise + log10 weight -------------------------------------

amr_sum_biomet <- amr_long_biomet %>%
  # rename BMD -> Antibiotic here so it propagates correctly
  mutate(Treatment = fct_recode(Treatment, "Antibiotic" = "BMD"),
         Treatment = as.character(Treatment)) %>%
  group_by(Treatment, Class) %>%
  summarise(
    total_abund = sum(Abundance),
    .groups = "drop"
  ) %>%
  mutate(weight = log10(total_abund + 1))

## ---- 3. Wide matrix (Treatment x AMR Class) ---------------------------

amr_mat_biomet <- amr_sum_biomet %>%
  select(Treatment, Class, weight) %>%
  pivot_wider(
    names_from  = Class,
    values_from = weight,
    values_fill = 0
  ) %>%
  column_to_rownames("Treatment") %>%
  as.matrix()

dim(amr_mat_biomet)
amr_mat_biomet[, 1:min(5, ncol(amr_mat_biomet))]

## ---- 4. Colors for treatments + AMR classes ---------------------------

amr_classes_biomet <- colnames(amr_mat_biomet)

amr_class_colors_biomet <- setNames(
  grDevices::colorRampPalette(c("#6c757d", "#adb5bd", "#ced4da"))(length(amr_classes_biomet)),
  amr_classes_biomet
)

grid_colors_biomet <- c(treat_colors, amr_class_colors_biomet)

## ---- 5. Long "links" table for chordDiagram ---------------------------

links_biomet <- amr_mat_biomet %>%
  as.data.frame() %>%
  rownames_to_column("Treatment") %>%
  pivot_longer(
    cols = -Treatment,
    names_to  = "Class",
    values_to = "weight"
  ) %>%
  filter(weight > 0)

head(links_biomet)

## ---- 6. Chord diagram: Treatments <-> Biocides & Metals ---------------

circos.clear()

chordDiagram(
  x = links_biomet,
  grid.col         = grid_colors_biomet,
  transparency     = 0.25,
  annotationTrack  = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    circos.text(
      x      = CELL_META$xcenter,
      y      = CELL_META$ylim[1] - mm_y(4),
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj    = c(0, 0.5),
      cex    = 0.5   # make labels small
    )
  },
  bg.border = NA
)

title("Chord diagram: Treatment \u2194 AMR biocide & metal classes (log10 abundance)")

## bring age into the plot 
library(dplyr)
library(tidyr)
library(tibble)
library(forcats)
library(circlize)

# Assumes:
# - ps_amr already defined
# - amr_long already created like you had:
#   amr_long <- psmelt(ps_amr) %>%
#       filter(Abundance > 0) %>%
#       filter(!is.na(Treatment), !is.na(Class))

BROADGROUP_COL <- "Broadclass"
BIO_MET_LABELS <- c("Biocides", "Metals")   # adjust if labels differ in your data

amr_long <- phyloseq::psmelt(ps_amr) %>%
  filter(Abundance > 0) %>%
  filter(!is.na(Treatment), !is.na(Class)) %>%
  mutate(
    Class      = gsub("_resistance", "", Class, ignore.case = TRUE) %>% trimws(),
    Broadclass = gsub("_resistance", "", Broadclass, ignore.case = TRUE) %>% trimws()
  )

# 1) Subset to biocides + metals, keep Days + Treatment
amr_long_biomet <- amr_long %>%
  filter(
    !is.na(.data[[BROADGROUP_COL]]),
    .data[[BROADGROUP_COL]] %in% BIO_MET_LABELS,
    !is.na(Class),
    !is.na(Age_Days)
  ) %>%
  # recode BMD -> Antibiotic
  mutate(
    Treatment = fct_recode(Treatment, "Antibiotic" = "BMD"),
    Treatment = as.character(Treatment),
    Days      = as.character(Age_Days)  # ensure we can paste it
  )

# 2) Summarise abundance by Treatment × Day × Class
amr_sum_biomet_age <- amr_long_biomet %>%
  group_by(Treatment, Age_Days, Class) %>%
  summarise(
    total_abund = sum(Abundance),
    .groups     = "drop"
  ) %>%
  mutate(
    weight      = log10(total_abund + 1),
    Trt_Day     = paste(Treatment, paste0("D", Age_Days), sep = "_")
    # e.g. "Basal Diet_D38"
  )

# 3) Wide matrix: rows = Treatment_Day, cols = AMR Class
amr_mat_biomet_age <- amr_sum_biomet_age %>%
  select(Trt_Day, Class, weight) %>%
  pivot_wider(
    names_from  = Class,
    values_from = weight,
    values_fill = 0
  ) %>%
  column_to_rownames("Trt_Day") %>%
  as.matrix()

dim(amr_mat_biomet_age)
amr_mat_biomet_age[, 1:min(5, ncol(amr_mat_biomet_age))]

# 4) Colors for Treatment×Day sectors and AMR classes ------------------

# Base treatment colors (as you had)
treat_colors <- c(
  "Basal Diet"     = "#eae69e",
  "Probiotic"      = "#bfdb81",
  "Essential oils" = "#83a561",
  "Antibiotic"     = "#48723e"
)

# Get all unique Trt_Day labels and their base Treatment
trt_day_labels <- rownames(amr_mat_biomet_age)
trt_day_treat  <- gsub("_D.*$", "", trt_day_labels)  # pull everything before "_D"

# Map each Trt_Day to its treatment color
trt_day_colors <- setNames(
  treat_colors[trt_day_treat],
  trt_day_labels
)

# AMR class colors (greyscale palette)
amr_classes_biomet <- colnames(amr_mat_biomet_age)

amr_class_colors_biomet <- setNames(
  grDevices::colorRampPalette(c("#6c757d", "#adb5bd", "#ced4da"))(length(amr_classes_biomet)),
  amr_classes_biomet
)

# Combine into one color vector for chordDiagram
grid_colors_biomet_age <- c(trt_day_colors, amr_class_colors_biomet)

# 5) Long links table ---------------------------------------------------

links_biomet_age <- amr_mat_biomet_age %>%
  as.data.frame() %>%
  rownames_to_column("Trt_Day") %>%
  pivot_longer(
    cols      = -Trt_Day,
    names_to  = "Class",
    values_to = "weight"
  ) %>%
  filter(weight > 0)

head(links_biomet_age)

# 6) Chord diagram: Treatment×Day <-> Biocide & Metal classes ----------

circos.clear()

chordDiagram(
  x               = links_biomet_age,
  grid.col        = grid_colors_biomet_age,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector_name <- CELL_META$sector.index
    circos.text(
      x      = CELL_META$xcenter,
      y      = CELL_META$ylim[1] - mm_y(4),
      labels = sector_name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj    = c(0, 0.5),
      cex    = 0.55  # slightly smaller, since more sectors
    )
  },
  bg.border = NA
)

title("Chord diagram: Treatment × Age \u2194 AMR biocide & metal classes (log10 abundance)")

