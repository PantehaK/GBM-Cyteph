## ============================================================================
## 01_repertoire_diversity_bootstrap.R
##
## Deep (bulk, e.g. Adaptive ImmunoSEQ-style) TCRB repertoire sequencing
## analysis. Because sequencing/template depth varies between samples,
## diversity metrics are compared at a common rarefied depth using
## multinomial resampling (bootstrap), rather than computed directly on raw
## template counts.
##
## For every sample this computes, at a fixed rarefaction depth:
##   - clonotype richness (number of unique rearrangements observed)
##   - Shannon entropy
##   - Pielou's evenness and clonality (1 - evenness)
## each with a bootstrap mean and 95% percentile CI.
##
## Input : data/bulk_tcrseq/*.tsv
##           One tab-separated file per sample, in the common Adaptive
##           Biotechnologies immunoSEQ export format (columns include at
##           least: amino_acid, frame_type, templates, rearrangement).
##           Filenames are expected to encode sample/patient/visit, e.g.
##           "<patient>_<visit>.tsv" -- adjust `parse_sample_name()` for your
##           own naming convention.
## Output: output/csv/repertoire_diversity_bootstrap_<project_tag>.csv
## ============================================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

source("../../00_config.R")

bulk_tcr_dir <- file.path(data_dir, "bulk_tcrseq")

# ---- Rarefaction / bootstrap settings ---------------------------------------
rarefaction_depth <- 50000   # set <= the lowest total_templates across samples
n_bootstrap_iter   <- 10000

# ---- 1. Load and clean per-sample repertoire files --------------------------
tcr_files <- list.files(bulk_tcr_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(tcr_files) == 0) {
  stop("No .tsv repertoire files found under ", bulk_tcr_dir)
}

# Adjust to match your own sample-naming convention, e.g. "<patient>_<visit>.tsv"
parse_sample_name <- function(filename) {
  base <- str_remove(basename(filename), "\\.tsv$")
  tibble(
    sample  = base,
    patient = str_extract(base, "^[^_]+"),
    visit   = str_extract(base, "(?<=_).+$")
  )
}

tcr_df <- map_dfr(tcr_files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    bind_cols(parse_sample_name(f)[rep(1, nrow(.)), ])
})

tcr_df_clean <- tcr_df %>%
  filter(!is.na(amino_acid), frame_type == "In")  # productive rearrangements only

depth_summary <- tcr_df_clean %>%
  group_by(sample, patient, visit) %>%
  summarise(total_templates = sum(templates, na.rm = TRUE), n_clonotypes = n(), .groups = "drop")

min_depth <- min(depth_summary$total_templates)
message("Lowest total template count across samples: ", min_depth,
        " (rarefaction_depth is currently set to ", rarefaction_depth, ")")
if (rarefaction_depth > min_depth) {
  warning("rarefaction_depth exceeds the lowest-depth sample; lower it to ",
          min_depth, " (or exclude that sample) so all samples are rarefied consistently.")
}

# ---- 2. Bootstrap richness / Shannon entropy / clonality --------------------
bootstrap_diversity <- function(df, depth = rarefaction_depth, n_iter = n_bootstrap_iter) {
  probs <- df$templates / sum(df$templates)

  boot_df <- replicate(n_iter, {
    sampled <- as.numeric(rmultinom(1, size = depth, prob = probs))
    sampled_nonzero <- sampled[sampled > 0]
    sampled_probs <- sampled_nonzero / sum(sampled_nonzero)

    S <- length(sampled_nonzero)                      # richness
    H <- -sum(sampled_probs * log(sampled_probs))      # Shannon entropy
    J <- H / log(S)                                    # Pielou evenness

    c(richness = S, shannon = H, clonality = 1 - J)
  }) %>% t() %>% as.data.frame()

  tibble(
    mean_richness = mean(boot_df$richness), richness_lower_ci = quantile(boot_df$richness, 0.025),
    richness_upper_ci = quantile(boot_df$richness, 0.975),
    mean_shannon = mean(boot_df$shannon), shannon_lower_ci = quantile(boot_df$shannon, 0.025),
    shannon_upper_ci = quantile(boot_df$shannon, 0.975),
    mean_clonality = mean(boot_df$clonality), clonality_lower_ci = quantile(boot_df$clonality, 0.025),
    clonality_upper_ci = quantile(boot_df$clonality, 0.975)
  )
}

diversity_results <- tcr_df_clean %>%
  group_by(sample, patient, visit) %>%
  group_split() %>%
  map_df(~ bind_cols(.x[1, c("sample", "patient", "visit")], bootstrap_diversity(.x)))

write.csv(diversity_results,
          file.path(csv_dir, paste0("repertoire_diversity_bootstrap_", project_tag, ".csv")),
          row.names = FALSE)
message("Wrote bootstrap diversity metrics for ", nrow(diversity_results), " sample(s).")
