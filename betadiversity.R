## Beta diversity 

#load packages 
#generally useful ones 
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)

#ones we specifically need 
library(phyloseq)
library(ggplot2)
library(vegan) # version 2.6-4
library(decontam)
library(microViz)
library(dplyr)
library(pairwiseAdonis)

#Load data
ps <- readRDS("ps/decontam-ps.rds")
ps

sample_data(ps)
#principal coordinate analysis

#colors
trts <- c("Basal Diet" = "#996666", 
          "Probiotic" = "#24868EFF",
          "Essential oils" = "#35B779FF",
          "BMD" = "#C7E020FF")

AGE <- c("1" ="#996699", 
         "10" = "#006666", 
         "21" = "#333399")

## Beta Diversity
ps <- tax_fix(ps)
ps <- phyloseq_validate(ps, remove_undetected = TRUE)

ps %>% 
  tax_transform(trans = "clr", rank = "Gene")

ps %>% 
  tax_transform(trans = "identity", rank = "unique") %>% 
  dist_calc("aitchison") # bray curtis distance

sample_data(ps)

# Plot a PCA
# Principal Components Analysis is an unconstrained method that does not use a distance matrix.
# Each point is a sample, and samples that appear closer together are typically 
# more similar to each other than samples which are further apart.


ps %>% tax_fix() %>% 
  tax_transform("clr", rank = "Gene") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Age_Days", shape = "Treatment", plot_taxa = 1:15, taxa_label_size = 7,tax_lab_style = tax_lab_style(size = 4, fontface = "bold.italic"), size = 6, alpha = 1.0) +
  scale_color_manual(values = c("#996699", "#006666", "#333399")) +
  scale_shape_manual(values = c(16, 15, 23, 17)) +
theme_classic() +
  labs(caption = "") +
  theme(axis.text= element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(hjust = 0.0)) +
  guides(fill = "none")

ggsave("plots/betadiversity_plot1.pdf")


# roughly estimate which samples will contain more of that taxon 
ps %>% tax_fix() %>% 
  tax_transform("clr", rank = "Gene") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Age_Days", shape = "Treatment", plot_taxa = 1:5, size = 4) +
  scale_color_manual(values = c("#996699", "#006666", "#333399")) +
  scale_shape_manual(values = c(16, 15, 23, 17)) +
  theme_classic()

ggsave(ggsave("plots/plot_PCA_age1-10-21.pdf", width = 25, height = 20, units = "cm"))

## PCoA
# Aitchison distance

ps %>% tax_fix() %>% 
  tax_transform("identity", rank = "Class") %>% # don't transform!
  dist_calc("aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment", size = 3) +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic"))


ggsave("plots/MDS_plot24.pdf")

#PERMANOVA overall trt

# calculate distances
aitchison_dists <- ps %>%
  tax_filter(min_prevalence = 0.1) %>%
  tax_transform("identity", rank = "Gene") %>%
  dist_calc("aitchison")

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Treatment"
  )

# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()

#Df  SumOfSqs         R2         F Pr(>F)
#Treatment  3   8748.54 0.02316144 0.7429326   0.95
#Residual  94 368971.44 0.97683856        NA     NA
#Total     97 377719.98 1.00000000        NA     NA

# view the info stored about the distance calculation
info_get(aitchison_perm)


#PERMANOVA overall Age (Days)
# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Gene") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment, permutations = 10000) 
# p = 0.99

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Age_Days, permutations = 10000)
# p = 9.999e-05 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Age_Days)

# calculate distances
aitchison_dists <- ps %>%
  tax_filter(min_prevalence = 0.1) %>%
  tax_transform("identity", rank = "Gene") %>%
  dist_calc("aitchison")

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Age_Days"
  )

# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()

pairwise.adonis(aitchison_dists, phyloseq::sample_data(ps)$Age_Days)


#Df  SumOfSqs        R2        F Pr(>F)
#Age_Days  2  91632.05 0.2425925 15.21393   0.01
#Residual 95 286087.93 0.7574075       NA     NA
#Total    97 377719.98 1.0000000       NA     NA

# view the info stored about the distance calculation
info_get(aitchison_perm)

#Centroids 
#subset by age and plot centroids by tRt
#d1
ps <- subset_samples(ps, Age_Days == "1")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Gene") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids


#Calculating centroids - June 2025
library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
  ) + 
  theme_classic()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function


#Centroids 
#subset by age and plot centroids by tRt
#d10
#Load data
ps <- readRDS("ps/decontam-ps.rds")
ps

ps <- subset_samples(ps, Age_Days == "10")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Gene") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids


#Calculating centroids - June 2025
library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
  ) + 
  theme_classic()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function



#Centroids 
#subset by age and plot centroids by tRt
#d21
#Load data
ps <- readRDS("ps/decontam-ps.rds")
ps

ps <- subset_samples(ps, Age_Days == "21")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Gene") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids


#Calculating centroids - June 2025
library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")

#Age centroids and distance lines 
plot_df_age <- left_join(plot_df, centroids_age, by = "Age_Days")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
  ) + 
  theme_classic()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function


#MDS plots 
# =========================
# NMDS + Hulls + Trajectories + PERMANOVA (updated)
# =========================

# --- Packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)   # for unit()

# --- Load phyloseq object
ps <- readRDS("ps/decontam-ps.RDS")  # <- keep your exact file name

# --- Colors
trts <- c(
  "Basal Diet"     = "#eae69e",
  "Probiotic"      = "#bfdb81",
  "Essential oils" = "#83a561",
  "Antibiotic"     = "#48723e"
)
ages <- c("1" = "#996699", "10" = "#006666", "21" = "#333399")

# =========================
# 0) Fix Treatment naming & factor levels (BMD -> Antibiotic)
# =========================
sd <- as.data.frame(sample_data(ps))
sd$Treatment <- trimws(sd$Treatment)
sd$Treatment[sd$Treatment == "BMD"] <- "Antibiotic"
sd$Treatment <- factor(sd$Treatment,
                       levels = c("Basal Diet","Probiotic","Essential oils","Antibiotic")
)
sample_data(ps)$Treatment <- sd$Treatment

# Rebuild Group (Treatment x Age) AFTER fixing Treatment names
sample_data(ps)$Group <- with(sample_data(ps),
                              interaction(Treatment, Age_Days, sep = " / ")
)

# =========================
# 1) CLR + Aitchison distance
# =========================
ps_clr <- ps %>%
  tax_fix() %>%
  tax_transform("clr", rank = "Gene") %>%
  ps_get()

aitch <- phyloseq::distance(ps_clr, method = "euclidean")  # Aitchison

# =========================
# 2) NMDS
# =========================
set.seed(123)
nmds <- vegan::metaMDS(as.matrix(aitch), k = 2, trymax = 200)

scores_df <- as.data.frame(nmds$points) %>%
  rownames_to_column("SampleID") %>%
  rename(MDS1 = MDS1, MDS2 = MDS2)

meta <- as(sample_data(ps_clr), "data.frame") %>%
  rownames_to_column("SampleID")

plot_df <- left_join(scores_df, meta, by = "SampleID")

# =========================
# 3) SAFE hulls (skip groups with <3 samples)
# =========================
safe_hulls <- function(df, grp_col) {
  df %>%
    group_by(.data[[grp_col]]) %>%
    filter(n() >= 3) %>%
    slice(chull(MDS1, MDS2)) %>%
    ungroup()
}

# -------------------------
# 3A) Treatment NMDS (Group hulls, colored by Treatment)
# -------------------------
hulls_trt <- safe_hulls(plot_df, "Group")

p_trt <- ggplot(plot_df, aes(MDS1, MDS2)) +
  geom_polygon(data = hulls_trt,
               aes(group = Group, color = Treatment, fill = Treatment),
               alpha = 0.30, linewidth = 0.6) +
  geom_point(aes(color = Treatment), size = 4) +
  scale_color_manual(values = trts,
                     breaks = names(trts),
                     labels = names(trts),
                     name   = "Treatment") +
  scale_fill_manual(values = trts,
                    breaks = names(trts),
                    labels = names(trts),
                    name   = "Treatment") +
  theme_classic() + labs(x="MDS1", y="MDS2",
                         title = "") +
  theme(
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.line    = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
p_trt


## no hulls MDS 
plot_df$Treatment <- factor(plot_df$Treatment, levels = names(trts))

p_mds_by_trt2 <- ggplot(plot_df, aes(MDS1, MDS2)) +
  geom_point(aes(color = Treatment), size = 3.5, alpha = 1.0) +
  stat_ellipse(aes(color = Treatment), linewidth = 1, linetype = 2) +  # 👈 add ellipses
  scale_color_manual(values = trts,
                     breaks = names(trts),
                     labels = names(trts)) +
  coord_equal() +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2") +
  theme(
    legend.position = "right",  # keep legend since it's useful for ellipses
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    strip.text   = element_text(size = 16, face = "bold"),
    axis.line    = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )

p_mds_by_trt2

# -------------------------
# 3B) Age NMDS (hulls & colors by Age)
# -------------------------
hulls_age <- plot_df %>%
  group_by(Age_Days) %>%
  filter(n() >= 3) %>%
  slice(chull(MDS1, MDS2)) %>%
  ungroup()

p_age <- ggplot(plot_df, aes(MDS1, MDS2)) +
  geom_polygon(data = hulls_age,
               aes(group = Age_Days, color = Age_Days, fill = Age_Days),
               alpha = 0.25, linewidth = 1.0) +
  geom_point(aes(color = Age_Days), size = 4) +
  scale_color_manual(values = ages, name = "Age (days)") +
  scale_fill_manual(values = ages,  name = "Age (days)") +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2",
       title = "") +
  theme(
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.line    = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )

p_age

# =========================
# 4) Centroid trajectories (Treatment across ages)
# =========================
centroids <- plot_df %>%
  group_by(Treatment, Age_Days) %>%
  summarise(cMDS1 = mean(MDS1), cMDS2 = mean(MDS2), .groups = "drop")

centroids$Age_Days <- factor(centroids$Age_Days, levels = c("1","10","21"))
centroids <- centroids %>% arrange(Treatment, Age_Days)

p_traj <- ggplot() +
  geom_point(data = plot_df, aes(MDS1, MDS2, color = Treatment), alpha = 0.35, size = 3) +
  geom_point(data = centroids, aes(cMDS1, cMDS2, color = Treatment), size = 6) +
  geom_path(data = centroids,
            aes(cMDS1, cMDS2, color = Treatment, group = Treatment),
            linewidth = 1.6, arrow = arrow(length = unit(0.25, "cm")),
            na.rm = TRUE) +
  geom_text(data = centroids,
            aes(cMDS1, cMDS2, label = Age_Days),
            size = 3, fontface = "bold", vjust = -1.1, color = "black") +
  scale_color_manual(values = trts, name = "Treatment") +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2",
       title = "",
       subtitle = "") +
  facet_wrap(~Treatment) + 
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )
p_traj

# =========================
# 5) PERMANOVA
# =========================
# 5A) Age only
meta_all <- as(sample_data(ps_clr), "data.frame")
set.seed(123)
adonis_age <- adonis2(aitch ~ Age_Days, data = meta_all, permutations = 999)

# 5B) Treatment only
set.seed(123)
adonis_trt <- adonis2(aitch ~ Treatment, data = meta_all, permutations = 999)

# 5C) Group = Treatment × Age
# (this matches your significant interaction model)
set.seed(123)
meta_all$Group <- interaction(meta_all$Treatment, meta_all$Age_Days, sep = " / ")
adonis_group <- adonis2(aitch ~ Group, data = meta_all, permutations = 999)

# 5D) Treatment within each age (your per-age tests)
per_age_results <- list()
for (age in c("1","10","21")) {
  ps_age <- subset_samples(ps_clr, Age_Days == age)
  dist_age <- phyloseq::distance(ps_age, method = "euclidean")
  meta_age <- as(sample_data(ps_age), "data.frame")
  # order alignment
  sids <- attr(dist_age, "Labels")
  meta_age <- meta_age[sids, , drop = FALSE]
  # ensure factor
  meta_age$Treatment <- droplevels(factor(meta_age$Treatment))
  set.seed(123)
  per_age_results[[paste0("Day_", age)]] <-
    adonis2(dist_age ~ Treatment, data = meta_age, permutations = 999)
}

# =========================
# 6) Print key outputs
# =========================
print(nmds$stress)
print(adonis_age)
print(adonis_trt)
print(adonis_group)
print(per_age_results)

# =========================
# 7) Save plots (optional)
# =========================
# dir.create("plots", showWarnings = FALSE, recursive = TRUE)
# ggsave("plots/NMDS_byTreatment_hulls.pdf", p_trt,  width = 8, height = 6)
# ggsave("plots/NMDS_byAge_hulls.pdf",       p_age,  width = 8, height = 6)
# ggsave("plots/NMDS_centroid_trajectories.pdf", p_traj, width = 8, height = 6)

# =========================
# 8) Quick sanity checks (optional)
# =========================
# Verify treatment counts by age (ensure Antibiotic is present where expected)
with(as.data.frame(sample_data(ps)), addmargins(table(Age_Days, Treatment)))
unique(plot_df$Treatment)

#interaction - Treatment*Age + Age
ps <- readRDS("ps/decontam-ps.rds")
ps

#PERMANOVA overall trt and age 
# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix() %>% 
  tax_transform(trans = "clr", rank = "Gene") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment

#interaction
sd <- as(phyloseq::sample_data(beta1), "data.frame")
colnames(sd)
sd$Age_Days <- factor(sd$Age_Days)
vegan::adonis2(
  psdist ~ Treatment * Age_Days + Age_Days,
  data = sd,
  permutations = 10000
)

## What are the treatments?


#Main effect of Days
vegan::adonis2(psdist ~ Age_Days, data = sd, permutations = 10000)

#Main effect of Treatment
vegan::adonis2(psdist ~ Treatment, data = sd, permutations = 10000)

#“Are there differences per treatment per age?”
vegan::adonis2(psdist ~ Treatment:Age_Days, data = sd, permutations = 10000)


## pairwise comparisons treatment per day 
library(purrr)
library(dplyr)
library(vegan)

# distance matrix
dist_mat <- as.matrix(psdist)

# all days
days_list <- sort(unique(sd$Age_Days))

pairwise_day_treat <- map_dfr(days_list, function(d) {
  
  # subset metadata to this day
  sub_sd_day <- sd %>% filter(Age_Days == d)
  
  # which treatments are present at this day?
  trts <- unique(sub_sd_day$Treatment)
  
  # if fewer than 2 treatments, skip this day
  if (length(trts) < 2) return(NULL)
  
  # all pairwise combinations of treatments
  trt_pairs <- combn(trts, 2, simplify = FALSE)
  
  # loop over treatment pairs within this day
  map_dfr(trt_pairs, function(tr) {
    
    sub_sd_pair <- sub_sd_day %>%
      filter(Treatment %in% tr)
    
    # make sure we really have BOTH treatments with samples
    if (n_distinct(sub_sd_pair$Treatment) < 2) return(NULL)
    if (min(table(sub_sd_pair$Treatment)) < 1) return(NULL)
    
    # subset distance matrix to these samples
    sub_dist_pair <- dist_mat[rownames(sub_sd_pair), rownames(sub_sd_pair)]
    
    # run PERMANOVA for this day and this pair
    res <- vegan::adonis2(
      sub_dist_pair ~ Treatment,
      data = sub_sd_pair,
      permutations = 10000
    )
    
    # return a tidy row
    tibble(
      Day    = d,
      Treat1 = tr[1],
      Treat2 = tr[2],
      F      = res$F[1],
      R2     = res$R2[1],
      p      = res$`Pr(>F)`[1]
    )
  })
})

pairwise_day_treat

pairwise_day_treat %>% 
  filter(p < 0.05)

library(dplyr)
library(ggplot2)
library(tidyr)

# Create pair label
pairwise_plot_df <- pairwise_day_treat %>%
  mutate(
    Pair = paste(Treat1, "vs", Treat2),
    neglog10p = -log10(p)
  )

# Plot
ggplot(pairwise_plot_df, aes(x = Day, y = neglog10p, color = Pair)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 1) +
  labs(
    title = "-log10(p) Pairwise PERMANOVA by Day",
    x = "Day",
    y = expression(-log[10](p)),
    color = "Treatment Pair"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )


## pairwise per day 
dist_mat <- as.matrix(psdist)
days_list <- unique(sd$Age_Days)

pairwise_by_day <- map_df(days_list, function(d) {
  
  sub_sd <- sd %>% filter(Age_Days == d)
  
  sub_dist <- dist_mat[rownames(sub_sd), rownames(sub_sd)]
  
  res <- vegan::adonis2(
    sub_dist ~ Treatment,
    data = sub_sd,
    permutations = 10000
  )
  
  data.frame(
    Day = d,
    F  = res$F[1],
    R2 = res$R2[1],
    p  = res$`Pr(>F)`[1]
  )
})

pairwise_by_day

#pairwise_by_day
#Day         F         R2         p
#1  10 1.1263788 0.10436154 0.1288871
#2   1 0.7394759 0.07106148 0.9652035
#3  21 0.9969240 0.09650524 0.4672533

## Treatment 
treat_list <- unique(sd$Treatment)

pairwise_by_treat <- map_df(treat_list, function(t) {
  
  sub_sd <- sd %>% filter(Treatment == t)
  
  sub_dist <- dist_mat[rownames(sub_sd), rownames(sub_sd)]
  
  res <- vegan::adonis2(
    sub_dist ~ Age_Days,
    data = sub_sd,
    permutations = 10000
  )
  
  data.frame(
    Treatment = t,
    F  = res$F[1],
    R2 = res$R2[1],
    p  = res$`Pr(>F)`[1]
  )
})

pairwise_by_treat

#pairwise_by_treat
#Treatment        F        R2         p
#1            BMD 5.945921 0.3615438 9.999e-05
#2     Basal Diet 5.001577 0.3030969 9.999e-05
#3      Probiotic 8.895085 0.4586257 9.999e-05
#4 Essential oils 4.933655 0.3196686 9.999e-05

#In other words:
#Treatments do not differ at Day 1
#Treatments do not differ at Day 10
#Treatments do not differ at Day 21
#But the Feed additives changed over time
#This produces the significant interaction

library(phyloseq)
library(ggplot2)
library(dplyr)

# ---- PCoA coordinates ----
ord <- ordinate(beta1, method = "PCoA", distance = psdist)
df_ord <- as.data.frame(ord$vectors)
df_ord$SampleID <- rownames(df_ord)

# Merge sample metadata
df_ord <- df_ord %>%
  left_join(sd %>% mutate(SampleID = rownames(sd)), by = "SampleID")

# ---- Compute centroids per Treatment × Day ----
centroids <- df_ord %>%
  group_by(Treatment, Age_Days) %>%
  summarise(
    PC1 = mean(Axis.1),
    PC2 = mean(Axis.2),
    PC3 = mean(Axis.3),
    .groups = "drop"
  ) %>%
  arrange(Treatment, Age_Days)


# ---- Plot ----
p12 <- ggplot() +
  geom_point(data = df_ord,
             aes(Axis.1, Axis.2, color = Treatment, shape = Age_Days),
             size = 5, alpha = 0.6) +
  geom_line(data = centroids,
            aes(PC1, PC2, color = Treatment, group = Treatment),
            size = 1.6) +
  theme_bw(base_size = 16) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2",
       title = "Trajectory: PC1 vs PC2 - AMR") +
  scale_color_manual(values = c(
    "Baseline"       = "#2ca02c",
    "Basal Diet"     = "#eae69e",
    "Probiotic"      = "#bfdb81",
    "Essential oils" = "#83a561",
    "BMD"     = "#48723e"
  )) +
  theme(panel.grid = element_blank())

p12 


p13 <- ggplot() +
  geom_point(data = df_ord,
             aes(Axis.1, Axis.3, color = Treatment, shape = Age_Days),
             size = 5, alpha = 0.6) +
  geom_line(data = centroids,
            aes(PC1, PC3, color = Treatment, group = Treatment),
            size = 1.6) +
  theme_bw(base_size = 16) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 3",
       title = "Trajectory: PC1 vs PC3") +
  scale_color_manual(values = c(
    "Baseline"       = "#2ca02c",
    "Basal Diet"     = "#eae69e",
    "Probiotic"      = "#bfdb81",
    "Essential oils" = "#83a561",
    "BMD"     = "#48723e"
  )) +
  theme(panel.grid = element_blank())

p13

p23 <- ggplot() +
  geom_point(data = df_ord,
             aes(Axis.2, Axis.3, color = Treatment, shape = Age_Days),
             size = 5, alpha = 0.6) +
  geom_line(data = centroids,
            aes(PC2, PC3, color = Treatment, group = Treatment),
            size = 1.6) +
  theme_bw(base_size = 16) +
  labs(x = "PCoA Axis 2", y = "PCoA Axis 3",
       title = "Trajectory: PC2 vs PC3") +
  scale_color_manual(values = c(
    "Baseline"       = "#2ca02c",
    "Basal Diet"     = "#eae69e",
    "Probiotic"      = "#bfdb81",
    "Essential oils" = "#83a561",
    "BMD"     = "#48723e"
  )) +
  theme(panel.grid = element_blank())

p23

