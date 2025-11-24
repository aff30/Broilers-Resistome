library(dplyr)
library(stringr)
library(viridis)
library(pheatmap)
library(tibble)
library(phyloseq)
#install.packages("microbiome")
library(microbiome)
library(tibble)
library(dplyr)
#install.packages("OTUtable")
library(OTUtable)
library(viridis)
library(stringr)
library(readxl) # if your metadata is an excel sheet and not csv
library(pheatmap)
library(ggplot2)


# load phyloseq object
ps <- readRDS("ps/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate the taxa at the desired level (Broadclass, Class, Mechanism, or Gene)
psmech <- aggregate_taxa(psrel, level = "Class")

# export otu and tax tables to use for heatmap
otu <- as.data.frame(otu_table(psmech))
tax <- as.data.frame(tax_table(psmech))

# convert the rownames to a column to merge the tables
otu2 <- rownames_to_column(otu, "Class")

# merge the tax table and otu tables by common variable (desired level)
merge <- merge(otu2, tax, "Class")

# order the merged table so that the classes are sorted by Broadclass (for aesthetics of the plot so they can be grouped with separation between broadclasses)
m2 <- merge[order(merge$Broadclass),]

# remove unwanted columns, only select the samples and the class
df2 <- dplyr::select(m2, !c("Broadclass", "unique"))

# remove the rownames so that the class can be assigned to the rownames in the desired order
df2 <- rownames_to_column(df2, "random")

# remove the unwanted columns, only select the samples and the class
df2 <- dplyr::select(df2, !c("random"))

# assign classes to rownames
merge2 <- column_to_rownames(df2, "Class")

# calculate z scores of the OTU table we edited
z <- OTUtable::zscore(merge2)

# convert z scores to a matrix
z3 <- as.matrix(z)

# use the tax table to define the classes to broadclass
# select only the columns we want
merge3 <- dplyr::select(merge, c("Class", "Broadclass"))

# convert to data frame
merge3 <- as.data.frame(merge3)

# assign the class to rownames
merge3 <- column_to_rownames(merge3, "Class")

# rename the column "Broadclass" for aesthetic purposes on the heatmap
names(merge3)[names(merge3) == "Broadclass"] <- "Resistance Type"

#fix rownames on the matrix for aestheic purposes on the heatmap 
rownames(z3) <- str_replace_all(rownames(z3), "_", " ")#sub underscores with spaces
rownames(z3) <- str_replace(rownames(z3), "resistance", "")#drop resistance because that's a given
rownames(z3) <- str_replace_all(rownames(z3), "and", "+") #simplfy this to make shorter

# edit rownames in the broadclass info table to match z3 matrix
rownames(merge3) <- str_replace_all(rownames(merge3), "_", " ")#sub underscores with spaces
rownames(merge3) <- str_replace(rownames(merge3), "resistance", "")#drop resistance because that's a given
rownames(merge3) <- str_replace_all(rownames(merge3), "and", "+") #simplfy this to make shorter

# --- After you run clean_names() on z3 and merge3 ---
fix_caps <- function(x) {
  # make "betalactams", "beta lactams", "beta-lactams" => "Betalactams"
  stringr::str_replace_all(
    x,
    stringr::regex("\\bbeta[ -]?lactams?\\b", ignore_case = TRUE),
    "Betalactams"
  )
}

# apply to row labels used in the heatmap
rownames(z3)     <- fix_caps(rownames(z3))
rownames(merge3) <- fix_caps(rownames(merge3))

# (if you later rebuild row_split, do it after these replacements)
rt_vec   <- as.character(merge3[,"Resistance Type", drop = TRUE])
row_split <- factor(rt_vec, levels = intersect(names(rt_cols), unique(rt_vec)))


# import metadata to assign variables at the top of the heatmap
# read in metadata
met <- read.delim("/Users/aff30/OneDrive - The Pennsylvania State University/Documents/P.h.D_projects/Year_1/amr++_trainning/AMR-tutorial/app-output/AMR_metadata.txt", sep = "\t")

# convert the sample IDs to rownames
met <- column_to_rownames(met, "Sample_ID")

# select only the columns we want
met2 <- dplyr::select(met, c("Treatment", "Age_Days", "Phase"))

# rename the columns for aesthetic purposes on the heatmap
names(met2)[names(met2) == "Phase"] <- "Feeding Phase"
names(met2)[names(met2) == "Age_Days"] <- "Age"

# --- 1) Clean rownames consistently in both z3 and merge3 ---
clean_names <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace("resistance", "") %>%
    str_replace_all("\\band\\b", "+") %>%
    str_squish()
}
rownames(z3)     <- clean_names(rownames(z3))
rownames(merge3) <- clean_names(rownames(merge3))

# --- 2) Read metadata and prep annotations ---
met <- read.delim("/Users/aff30/OneDrive - The Pennsylvania State University/Documents/P.h.D_projects/Year_1/amr++_trainning/AMR-tutorial/app-output/AMR_metadata.txt", sep = "\t") |>
  column_to_rownames("Sample_ID") |>
  dplyr::select(Treatment, Age_Days, Phase) |>
  dplyr::rename(`Feeding Phase` = Phase, Age = Age_Days)

# Optional rename for legend consistency
met$Treatment <- dplyr::recode(met$Treatment, "BMD" = "Antibiotic")


# Make Age numeric when possible (for ordering), then keep a character copy for colors
suppressWarnings({ age_num <- as.integer(as.character(met$Age)) })
met$Age <- ifelse(!is.na(age_num), age_num, met$Age)

# --- 3) Align z3 (columns) and metadata (rows) ---
common <- intersect(colnames(z3), rownames(met))
z3     <- z3[, common, drop = FALSE]
met2   <- met[common, c("Age","Treatment","Feeding Phase"), drop = FALSE]

# Order by Age (numeric if possible)
order_idx <- if (is.numeric(met2$Age)) order(met2$Age) else order(as.character(met2$Age))
z3   <- z3[, order_idx, drop = FALSE]
met2 <- met2[order_idx, , drop = FALSE]

# Convert Age to character so manual colors match exactly
met2$Age <- as.character(met2$Age)

# --- 4) Row annotations: Resistance Type alignment + row gaps ---
merge3 <- merge3[rownames(z3), , drop = FALSE]
rt <- merge3[["Resistance Type"]]
row_order <- order(rt)
z3       <- z3[row_order, , drop = FALSE]
merge3   <- merge3[row_order, , drop = FALSE]
rt_sorted <- merge3[["Resistance Type"]]
gaps_row <- which(rt_sorted[-1] != rt_sorted[-length(rt_sorted)])

# --- 5) Annotation color palettes ---
## Age: your exact custom colors (assumes ages 1, 10, 21)
age_palette <- c(
  "1"  = "#996699",
  "10" = "#006666",
  "21" = "#333399"
)

## Treatment: fixed (subset to present)
treat_palette <- c(
  "Baseline"       = "#2ca02c",
  "Basal Diet"     = "#eae69e",
  "Probiotic"      = "#bfdb81",
  "Essential oils" = "#83a561",
  "Antibiotic"     = "#48723e"
)
treat_vals <- unique(as.character(met2$Treatment))
treat_cols <- treat_palette[intersect(names(treat_palette), treat_vals)]

## Feeding Phase: custom distinct colors (subset to present)
phase_palette <- c(
  "Starter"  = "#D35127FF",
  "Grower"   = "#272A59FF"
  # add "Finisher" = "#3B6FB6" if you ever have it
)
phase_vals <- unique(as.character(met2[["Feeding Phase"]]))
phase_cols <- phase_palette[intersect(names(phase_palette), phase_vals)]
if (length(phase_cols) < length(phase_vals)) {
  missing <- setdiff(phase_vals, names(phase_cols))
  auto    <- setNames(viridis::magma(length(missing)), missing)
  phase_cols <- c(phase_cols, auto)
}

## Resistance Type: fixed + fallback if needed
rt_palette <- c(
  "Biocides"       = "#552F7A",
  "Drugs"          = "#7C5F98",
  "Metals"         = "#B09FC1",
  "Multi-compound" = "#CABED6"
)
rt_vals <- unique(as.character(rt_sorted))
rt_cols <- rt_palette[intersect(names(rt_palette), rt_vals)]
if (any(!rt_vals %in% names(rt_palette))) {
  missing <- setdiff(rt_vals, names(rt_palette))
  rt_cols <- c(rt_cols, setNames(viridis::plasma(length(missing)), missing))
}

annotation_colors <- list(
  "Age"             = age_palette,
  "Treatment"       = treat_cols,
  "Resistance Type" = rt_cols
)

# --- 6) Draw the heatmap (dark-purple gradient) ---
dark_purple <- colorRampPalette(c("#F2F2F2FF", "#7C7189FF", "#FAE093FF", "#D04E59FF", "#BC8E7DFF", "#2F3D70FF"))(30)

heat <- pheatmap(
  z3,
  annotation_col       = subset(met2, select = c("Age","Treatment")),
  annotation_row       = subset(merge3, select = "Resistance Type"),
  cluster_rows         = FALSE,
  cluster_cols         = FALSE,
  color                = dark_purple,
  annotation_colors    = annotation_colors,
  gaps_row             = gaps_row,
  legend               = TRUE,
  annotation_legend    = TRUE,
  show_rownames        = TRUE,
  show_colnames        = FALSE,
  border_color         = NA,
  annotation_names_col = TRUE
)

heat

### legend in bottom 
# --- Libraries ---
library(dplyr)
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Assume you already have: z3 (matrix of z-scores) and merge3 (data.frame with "Resistance Type") ---

# ---- 1) Metadata (Age & Treatment only) ----
# Load metadata and align to z3 as you were doing
met <- read.delim("/Users/aff30/OneDrive - The Pennsylvania State University/Documents/P.h.D_projects/Year_1/amr++_trainning/AMR-tutorial/app-output/AMR_metadata.txt", sep = "\t") |>
  column_to_rownames("Sample_ID") |>
  dplyr::select(Treatment, Age_Days) |>
  dplyr::rename(Age = Age_Days)

# Optional rename for legend consistency
met$Treatment <- dplyr::recode(met$Treatment, "BMD" = "Antibiotic")

# numeric age for ordering; keep character for annotation coloring
suppressWarnings({ age_num <- as.integer(as.character(met$Age)) })
met$Age <- ifelse(!is.na(age_num), age_num, met$Age)

# Align columns (samples) between z3 and metadata
common <- intersect(colnames(z3), rownames(met))
z3     <- z3[, common, drop = FALSE]
met2   <- met[common, c("Age","Treatment"), drop = FALSE]

# Order columns by Age if numeric
order_idx <- if (is.numeric(met2$Age)) order(met2$Age) else order(as.character(met2$Age))
z3   <- z3[, order_idx, drop = FALSE]
met2 <- met2[order_idx, , drop = FALSE]

# Convert Age to character for exact color name matching
met2$Age       <- as.character(met2$Age)
met2$Treatment <- as.character(met2$Treatment)

# ---- 2) Row annotation (Resistance Type) prep ----
# merge3 must be in same row order as z3
merge3 <- merge3[rownames(z3), , drop = FALSE]
rt_vec <- as.character(merge3[,"Resistance Type", drop = TRUE])

# ---- 3) Colors ----
# Heatmap (dark purple)
dark_purple <- colorRampPalette(c("#212E52FF", "#444E7EFF", "#8087AAFF", "#B7ABBCFF","#F7F7F7FF", "#F9ECE8FF", "#FCC893FF", "#FEB424FF", "#FD8700FF"))(30)

dark_purple <- colorRampPalette(c("#2166AC", "#BFD3E6", "#FFFFBF", "#B2182B"))(30)


col_fun <- colorRamp2(seq(min(z3, na.rm = TRUE), max(z3, na.rm = TRUE), length.out = length(dark_purple)),
                      dark_purple)

# Age (your exact 3 colors)
age_cols <- c(
  "1"  = "#996699",
  "10" = "#006666",
  "21" = "#333399"
)

# Treatment
treat_cols <- c(
  "Baseline"       = "#2ca02c",
  "Basal Diet"     = "#eae69e",
  "Probiotic"      = "#bfdb81",
  "Essential oils" = "#83a561",
  "Antibiotic"     = "#48723e"
)

# Resistance Type
rt_cols <- c(
  "Biocides"       = "#C67B6FFF",
  "Drugs"          = "#DE9B71FF",
  "Metals"         = "#EFBC82FF",
  "Multi-compound" = "#FBDFA2FF"
)
# (auto-assign any unseen categories)
if (any(!unique(rt_vec) %in% names(rt_cols))) {
  missing <- setdiff(unique(rt_vec), names(rt_cols))
  extra   <- setNames(colorRampPalette(c("#A29BFE","#D6CCFF","#7E57C2","#B39DDB"))(length(missing)), missing)
  rt_cols <- c(rt_cols, extra)
}

# ---- 4) Annotations ----
top_anno <- HeatmapAnnotation(
  Age       = met2$Age,
  Treatment = met2$Treatment,
  col = list(
    Age       = age_cols,
    Treatment = treat_cols
  ),
  annotation_name_side = "right",
  annotation_name_gp   = gpar(fontsize = 10, fontface = "bold"),  # <-- makes Age & Treatment bold
  gp = gpar(col = NA)  # no borders
)

# Hide the sidebar labels for Resistance Type (keep legend)

left_anno <- rowAnnotation(
  "Resistance Type" = rt_vec, #resistance type here is a vector
  col = list("Resistance Type" = rt_cols),
  show_annotation_name = FALSE,         # <- hides "Biocides / Drugs / ..." text on the side
  annotation_name_gp   = gpar(fontsize = 0),
  annotation_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  gp = gpar(col = NA)
)

# Optional: split rows by Resistance Type (gives subtle group bands)
row_split <- factor(rt_vec, levels = intersect(names(rt_cols), unique(rt_vec))) 


# ---- 5) Heatmap ----
## with dendograms 
ht <- Heatmap(
  z3,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,                #turn OFF row dendrograms
  cluster_columns = TRUE,              #keep column dendrograms
  clustering_distance_columns = "pearson",
  clustering_method_columns = "average",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_split = row_split,
  row_title = NULL,
  top_annotation = top_anno,
  row_names_gp = gpar(fontsize = 8),
  column_dend_height = unit(12, "mm"), # keep only column dendrogram sizing
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

draw(left_anno + ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")


# subset drugs
# ---- SUBSET: keep only rows where Resistance Type == "Drugs" ----
drug_idx <- rt_vec == "Drugs"

z3_drug   <- z3[drug_idx, , drop = FALSE]     # subset matrix
rt_vec_drug <- rt_vec[drug_idx]              # subset resistance type vector
merge3_drug <- merge3[drug_idx, , drop = FALSE]

# (Optional) if you have a Class or AMR_Class column you want as row names:
# rownames(z3_drug) <- merge3_drug$Class   # only if needed/desired

# ---- 4) Annotations (unchanged for columns) ----
top_anno <- HeatmapAnnotation(
  Age       = met2$Age,
  Treatment = met2$Treatment,
  col = list(
    Age       = age_cols,
    Treatment = treat_cols
  ),
  annotation_name_side = "right",
  annotation_name_gp   = gpar(fontsize = 10, fontface = "bold"),
  gp = gpar(col = NA)
)

# Row annotation is now trivial (only Drugs, but we keep legend if you like)
left_anno <- rowAnnotation(
  "Resistance Type" = rt_vec_drug,
  col = list("Resistance Type" = rt_cols),
  show_annotation_name = FALSE,
  annotation_name_gp   = gpar(fontsize = 0),
  annotation_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  gp = gpar(col = NA)
)

# No need for row_split now (all are Drugs), but if you leave it:
# row_split <- factor(rt_vec_drug, levels = intersect(names(rt_cols), unique(rt_vec_drug)))
row_split <- NULL

ht <- Heatmap(
  z3_drug,                     # <-- use z3_drug here
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  clustering_distance_columns = "pearson",
  clustering_method_columns = "average",
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_split = row_split,       # NULL is fine
  row_title = NULL,
  top_annotation = top_anno,
  row_names_gp = gpar(fontsize = 10),
  column_dend_height = unit(12, "mm"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

draw(left_anno + ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right")

