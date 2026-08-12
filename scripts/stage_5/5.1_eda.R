# Pre-Modelling EDA: NirK / PCuAC / Both final species-level co-occurrence 
# dataset

# Purpose: 
# Dataset-level exploratory data analysis performed BEFORE the dataset is 
# split into the downstream modelling questions:

# Model A: among species containing NirK, can PCuAC co-occurrence presence be predicted? 
# Model A: among species containing PCuAC, can NirK co-occurrence presence be predicted? 

# Unit of analysis:
# One row = one species-level observation
# Species-level observations were produced by the upstream species-level 
# presence/absence scripts (4.0-4.2). Sequence-derived features were calculated 
# once per species form a single representative accession per protein 
# (Bio.SeqUtils.ProtParam.ProteinAnalysis).


# Scope of this script

# DOES
# - inspect structure, size and class composition of the dataset
# - asses missingness, duplication, and general data quality
# - visualize NirK and PCuAC feature distributions (all 25 features each:
#   length, molecular weight, isoelectric point, GRAVY, aromaticity, and
#   the 20 raw amino-acid composition fractions)
# - compare feature values across classes
# - examine correlation / redundancy among features

# DOES NOT
# - train Logistic Regression or Random Forest models
# - perform train/test splitting, cross-validation, or feature scaling
#   for modelling
# - perform model-specific feature selection / correlation filtering
# - calculate model coefficients, permutation importance, or any
#   model-derived quantity

# Ctrl + Shift + R for sections




# Setup -------------------------------------------------------------------
library(tidyverse)
set.seed(11)

input_file <- "final_feature_table.csv"
output_dir <- "out/eda"
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Identifier columns 
metadata_cols <- c(
  "species", "taxid", "organism_name",
  "nirK_present", "pcuac_present", "class",
  "nirK_rep_accession", "pcuac_rep_accession", "species_resolved"
)

# Biological class labels
expected_classes <- c("nirK_only", "pcuac_only", "both")

# The 20 standard amino acids
aa_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S",
             "T","V","W","Y")

# Physicochemical features
# one set per protein
nirk_core_features <- c("nirK_length", "nirK_molecular_weight",
                        "nirK_isoelectric_point", "nirK_gravy",
                        "nirK_aromaticity")
pcuac_core_features <- c("pcuac_length", "pcuac_molecular_weight",
                        "pcuac_isoelectric_point", "pcuac_gravy",
                        "pcuac_aromaticity")

# Amino acid composition fractions
nirk_aa_features <- paste0("nirK_aa_frac_", aa_list)
pcuac_aa_features <- paste0("pcuac_aa_frac_", aa_list)

# Full feature sets used for distributions & correlations
nirk_features <- c(nirk_core_features, nirk_aa_features)
pcuac_features <- c(pcuac_core_features, pcuac_aa_features)

# Labels for features in plot titles/axes
feature_labels <- c(
  length = "Length",
  molecular_weight = "MW",
  isoelectric_point = "pI",
  gravy = "GRAVY",
  aromaticity = "Aromaticity"
)

# Subset used for group-wise comparisons
nirk_group_features <- c(nirk_core_features, "nirK_aa_frac_C", "nirK_aa_frac_H", 
                         "nirK_aa_frac_M")
pcuac_group_features <- c(pcuac_core_features, "pcuac_aa_frac_C", "pcuac_aa_frac_H", 
                          "pcuac_aa_frac_M")

stopifnot(all(nirk_features %in% c(nirk_core_features, nirk_aa_features)))
stopifnot(all(pcuac_features %in% c(pcuac_core_features, pcuac_aa_features)))




# Dataset loading and structure -------------------------------------------
df <- read_csv(input_file, show_col_types = FALSE)

missing_expected <- setdiff(c(metadata_cols, nirk_features, pcuac_features), names(df)) 
if (length(missing_expected) > 0) {
  stop("ERROR: expected column(s) not found in ", input_file, ": ",
       paste(missing_expected, collapse = ", "))
}

col_types <- tibble(
  column = names(df),
  type = map_chr(df, ~ class(.x)[1])
)
write_csv(col_types, file.path(tables_dir, "column_types.csv"))

dataset_structure_txt <- file.path(tables_dir, "dataset_structure.txt")
sink(dataset_structure_txt)
cat("Rows:", nrow(df), "\n")
cat("Columns:", ncol(df), "\n\n")
cat("Column names and types:\n")
print(col_types, n = Inf)
cat("\nSummary statistics (numeric columns):\n")
print(summary(df %>% select(where(is.numeric))))
sink()

# Class categories present in the final dataset 
class_counts_all <- df %>% count(class, name = "n") %>% 
  mutate(percent = round(100 * n / sum(n), 2))
write_csv(class_counts_all, file.path(tables_dir, "class_counts_full_dataset.csv"))

# Check for unexpected class labels
unexpected_classes <- setdiff(unique(df$class), expected_classes)
if (length(unexpected_classes) > 0) {
  warning("Unexpected class label(s) found and retained: ",
          paste(unexpected_classes, collapse = ", "))
}




# Data-quality (missing values & duplicates) ------------------------------
# Missing values:
missing_by_column <- function(cols, block_label) {
  tibble(
    block = block_label,
    column = cols,
    n_missing = map_int(df[cols], ~ sum(is.na(.x))),
  ) %>%
    mutate(percent_missing = round(100 * n_missing / nrow(df), 2))
}

missing_metadata <- missing_by_column(metadata_cols, "metadata")
missing_nirk <- missing_by_column(nirk_features, "nirK_feature")
missing_pcuac <- missing_by_column(pcuac_features, "pcuac_feature")

missing_summary <- bind_rows(missing_metadata, missing_nirk, missing_pcuac)
write_csv(missing_summary, file.path(tables_dir, "missingness_by_column.csv"))

# Structural missingness is expected for absent proteins
missingness_by_class <- df %>%
  mutate(
    nirK_features_missing = if_any(all_of(nirk_features), is.na),
    pcuac_features_missing = if_any(all_of(pcuac_features), is.na)
  ) %>%
  count(class, nirK_features_missing, pcuac_features_missing)
write_csv(missingness_by_class, file.path(tables_dir, "missingness_by_class.csv"))

# Flag unexpected missing values inconsistent with class
unexpected_missingness <- df %>%
  mutate(
    nirK_features_missing = if_any(all_of(nirk_features), is.na),
    pcuac_features_missing = if_any(all_of(pcuac_features), is.na)
  ) %>%
  filter(
    (class == "both" & (nirK_features_missing | pcuac_features_missing)) |
      (class == "nirK_only"  & nirK_features_missing) |
      (class == "pcuac_only" & pcuac_features_missing)
  ) %>%
  select(species, class, nirK_features_missing, pcuac_features_missing)

write_csv(unexpected_missingness, file.path(tables_dir, "unexpected_missingness.csv"))

if (nrow(unexpected_missingness) > 0) {
  warning(nrow(unexpected_missingness),
          " species have missingness not explained by class-structural absence ",
          "(see unexpected_missingness.csv) and REQUIRES MANUAL REVIEW.")
}


# Duplicates
n_duplicated_species <- sum(duplicated(df$species))

duplicate_feature_vectors <- function(features) {
  complete_df <- df %>% filter(if_all(all_of(features), ~ !is.na(.x)))
  is_dup <- duplicated(complete_df %>% select(all_of(features))) |
    duplicated(complete_df %>% select(all_of(features)), fromLast = TRUE)
  complete_df %>% filter(is_dup) %>% select(species, class, all_of(features))
}

dup_nirk_vectors <- duplicate_feature_vectors(nirk_features)
dup_pcuac_vectors <- duplicate_feature_vectors(pcuac_features)

n_dup_nirk_rows <- nrow(dup_nirk_vectors)
n_dup_pcuac_rows <- nrow(dup_pcuac_vectors)

n_dup_full_rows <- sum(duplicated(df) | duplicated(df, fromLast = TRUE))

write_csv(dup_nirk_vectors, file.path(tables_dir, "duplicate_nirk_feature_vectors.csv"))
write_csv(dup_pcuac_vectors, file.path(tables_dir, "duplicate_pcuac_feature_vectors.csv"))

# In-class vs across-class breakdown for each duplicate feature-vector set
dup_class_breakdown <- function(dup_df, block_label) {
  if (nrow(dup_df) == 0) return(tibble(block = block_label, note = "no duplicates found"))
  dup_df %>% count(class, name = "n_duplicated_rows") %>% mutate(block = block_label)
}
dup_class_summary <- bind_rows(dup_class_breakdown(dup_nirk_vectors, "nirK_feature_vector"),
                               dup_class_breakdown(dup_pcuac_vectors, "pcuac_feature_vector")
)
write_csv(dup_class_summary, file.path(tables_dir, "duplicate_class_breakdown.csv"))
# Note: duplicates are reported, not removed




# Class balance -----------------------------------------------------------
ggplot(class_counts_all, aes(x = class, y = n, fill = class)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (",percent,"%)")), vjust = -0.3) +
  labs(title = "Species by class (full dataset)", x = "", y = "Number of species") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "class_counts_bar.png"), width = 6, height = 5)

# Model A and Model B populations
model_a_population <- df %>% filter(class %in% c("nirK_only", "both"))
model_b_population <- df %>% filter(class %in% c("pcuac_only", "both"))

model_a_counts <- model_a_population %>% count(class, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 2), population = "Model A (nirK_only vs both)")

model_b_counts <- model_b_population %>% count(class, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 2), population = "Model B (pcuac_only vs both)")

modelling_population_counts <- bind_rows(model_a_counts, model_b_counts)

write_csv(modelling_population_counts, file.path(tables_dir, "modelling_population_counts.csv"))

ggplot(modelling_population_counts, aes(x = class, y = n, fill = class)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", percent, "%)")), vjust = -0.3) +
  facet_wrap(~population, scales = "free_x") +
  labs(title = "Model A / Model B population sizes",
       x = "", y = "Number of species") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(plots_dir, "modelling_population_counts_bar.png"), width = 8, height = 5)




# Protein feature distributions (physicochemical features & AA composition)--------

# Physicochemical features
core_long <- function(features, prefix, protein_label) {
  df %>%
    select(species, class, all_of(features)) %>%
    filter(if_all(all_of(features), ~ !is.na(.x))) %>%
    pivot_longer(all_of(features), names_to = "feature", values_to = "value") %>%
    mutate(
      feature_label = feature_labels[sub(paste0("^", prefix, "_"), "", feature)],
      protein = protein_label
    )
}

nirk_core_long <- core_long(nirk_core_features, "nirK", "NirK")
pcuac_core_long <- core_long(pcuac_core_features, "pcuac", "PCuAC")

plot_core_distribution <- function(long_df, protein_label, filename) {
  ggplot(long_df, aes(x = value, fill = class)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~feature_label, scales = "free", ncol = 3) +
    labs(title = paste(protein_label, "core feature distributions, by class"),
         x = "", y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())
  ggsave(file.path(plots_dir, filename), width = 10, height = 7)
}

plot_core_distribution(nirk_core_long, "NirK", "nirk_core_feature_distributions.png")
plot_core_distribution(pcuac_core_long, "PCuAC", "pcuac_core_feature_distributions.png")

# Amino-acid composition
aa_long <- function(features, prefix, protein_label) {
  df %>%
    select(species, class, all_of(features)) %>%
    filter(if_all(all_of(features), ~ !is.na(.x))) %>%
    pivot_longer(all_of(features), names_to = "feature", values_to = "value") %>%
    mutate(
      amino_acid = sub(paste0("^", prefix, "_aa_frac_"), "", feature),
      protein = protein_label
    )
}

nirk_aa_long <- aa_long(nirk_aa_features, "nirK", "NirK")
pcuac_aa_long <- aa_long(pcuac_aa_features, "pcuac", "PCuAC")

plot_aa_composition <- function(long_df, protein_label, filename) {
  ggplot(long_df, aes(x = class, y = value, fill = class)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~amino_acid, scales = "free_y", ncol = 5) +
    labs(title = paste(protein_label, "amino-acid composition, by class"),
         x = "", y = "Fraction") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom", legend.title = element_blank())
  ggsave(file.path(plots_dir, filename), width = 12, height = 9)
}

plot_aa_composition(nirk_aa_long, "NirK", "nirk_aa_composition_by_class.png")
plot_aa_composition(pcuac_aa_long, "PCuAC", "pcuac_aa_composition_by_class.png")

low_variance_features <- function(features, protein_label) {
  df %>%
    select(all_of(features)) %>%
    summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "sd") %>%
    mutate(protein = protein_label) %>%
    arrange(sd)
}

variance_summary <- bind_rows(
  low_variance_features(nirk_features,  "NirK"),
  low_variance_features(pcuac_features, "PCuAC")
)
write_csv(variance_summary, file.path(tables_dir, "feature_variance_summary.csv"))
# Low variance features are sorted at the top for review before modelling.




# Group-wise comparisons --------------------------------------------------
group_box_plot <- function(population, features, prefix, group_levels, title, filename) {
  long_df <- population %>%
    mutate(group = factor(class, levels = group_levels)) %>%
    select(group, all_of(features)) %>%
    filter(if_all(all_of(features), ~ !is.na(.x))) %>%
    pivot_longer(-group, names_to = "feature", values_to = "value") %>%
    mutate(
      feature_label = if_else(
        feature %in% c(nirk_core_features, pcuac_core_features),
        feature_labels[sub(paste0("^", prefix, "_"), "", feature)],
        paste0("%", toupper(sub(paste0("^", prefix, "_aa_frac_"), "", feature)))
      )
    )
  
  ggplot(long_df, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    facet_wrap(~feature_label, scales = "free_y", ncol = 4) +
    labs(title = title, x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "none")
  ggsave(file.path(plots_dir, filename), width = 10, height = 7)
}

group_box_plot(
  population = model_a_population,
  features = nirk_group_features,
  prefix = "nirK",
  group_levels = c("nirK_only", "both"),
  title = "NirK features: nirK_only vs both",
  filename = "modelA_group_comparison_boxplot.png"
)

group_box_plot(
  population = model_b_population,
  features = pcuac_group_features,
  prefix = "pcuac",
  group_levels = c("pcuac_only", "both"),
  title = "PCuAC features: pcuac_only vs both",
  filename = "modelB_group_comparison_boxplot.png"
)




# Correlation and feature redundancy --------------------------------------
correlation_matrix_and_pairs <- function(features, protein_label, prefix,
                                         heatmap_file, pairs_threshold = 0.8) {
  cor_data <- df %>% select(all_of(features)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
  
  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column("feature1") %>%
    pivot_longer(-feature1, names_to = "feature2", values_to = "correlation")
  
  p <- ggplot(cor_long, aes(x =feature1, y = feature2, fill= correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = paste(protein_label, "feature correlations (n =", nrow(cor_data), 
                       "complete species)"), x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6))
  ggsave(heatmap_file, plot = p, width = 10, height = 9)
  
  high_pairs <- cor_long %>%
    filter(feature1 < feature2, abs(correlation) >= pairs_threshold) %>%
    arrange(desc(abs(correlation)))
  
  list(matrix = cor_mat, high_pairs = high_pairs)
}

nirk_cor <- correlation_matrix_and_pairs(
  nirk_features, "NirK", "nirk",
  file.path(plots_dir, "nirk_correlation_heatmap.png")
)
pcuac_cor <- correlation_matrix_and_pairs(
  pcuac_features, "PCuAC", "pcuac",
  file.path(plots_dir, "pcuac_correlation_heatmap.png")
)

write_csv(as.data.frame(nirk_cor$matrix) %>% rownames_to_column("feature"),
          file.path(tables_dir, "nirk_correlation_matrix.csv"))
write_csv(as.data.frame(pcuac_cor$matrix) %>% rownames_to_column("feature"),
          file.path(tables_dir, "pcuac_correlation_matrix.csv"))

high_corr_pairs <- bind_rows(
  nirk_cor$high_pairs %>% mutate(protein = "NirK"),
  pcuac_cor$high_pairs %>% mutate(protein = "PCuAC")
) %>% relocate(protein)

write_csv(high_corr_pairs, file.path(tables_dir, "highly_correlated_feature_pairs.csv"))




# Print entire summary (y/n) ----------------------------------------------
eda_summary <- tribble(
  ~metric, ~value,
  "Total species observations", as.character(nrow(df)),
  "Total dataset columns", as.character(ncol(df)),
  "NirK-only count", as.character(class_counts_all$n[class_counts_all$class == "nirK_only"]),
  "PCuAC-only count", as.character(class_counts_all$n[class_counts_all$class == "pcuac_only"]),
  "both (NirK + PCuAC) count", as.character(class_counts_all$n[class_counts_all$class == "both"]),
  "NirK-only (%)", paste0(class_counts_all$percent[class_counts_all$class == "nirK_only"], "%"),
  "PCuAC-only (%)", paste0(class_counts_all$percent[class_counts_all$class == "pcuac_only"], "%"),
  "both NirK + PCuAC (%)", paste0(class_counts_all$percent[class_counts_all$class == "both"], "%"),
  "Unexpected class labels present", as.character(length(unexpected_classes) > 0),
  "Maximum metadata missingness (%)", as.character(max(missing_metadata$percent_missing)),
  "Maximum NirK feature missingness (%)", as.character(max(missing_nirk$percent_missing)),
  "Maximum PCuAC feature missingness (%)", as.character(max(missing_pcuac$percent_missing)),
  "Species with unexpected missingness", as.character(nrow(unexpected_missingness)),
  "Observations with duplicated species identifiers", as.character(n_duplicated_species),
  "Observations with duplicated NirK feature vector", as.character(n_dup_nirk_rows),
  "Observations with duplicated PCuAC feature vector", as.character(n_dup_pcuac_rows),
  "Fully duplicated observations", as.character(n_dup_full_rows),
  "Number of NirK features", as.character(length(nirk_features)),
  "Number of PCuAC features", as.character(length(pcuac_features)),
  "Highly correlated feature pairs (|r| >= 0.8)", as.character(nrow(high_corr_pairs)),
  "Model A population size", as.character(nrow(model_a_population)),
  "Model B population size", as.character(nrow(model_b_population))
)
write_csv(eda_summary, file.path(tables_dir, "eda_summary.csv"))

if (readline("Would you like to print EDA script summary: y/n? ") == "y") {
  print(eda_summary, nrow(eda_summary))
}














