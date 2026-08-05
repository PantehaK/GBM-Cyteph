## ============================================================================
## 02_clonal_expansion_significance_testing.R
##
## Per-patient, per-clonotype exact binomial test for significant clonal
## expansion/contraction between two timepoints (e.g. pre- vs post-product
## infusion), analogous to a paired differential-abundance test. For each
## clonotype detected in either timepoint, the number of templates observed
## post-timepoint is compared to a binomial expectation based on each
## timepoint's total sequencing depth; p-values are Benjamini-Hochberg
## corrected within each patient (and optionally globally).
##
## Input : data/bulk_tcrseq/*.tsv  (see 01_repertoire_diversity_bootstrap.R
##           for expected format and `parse_sample_name()`)
## Output: output/csv/clonal_expansion_significant_<project_tag>.csv
##         output/csv/clonal_expansion_all_tested_<project_tag>.csv
## ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)

source("../../00_config.R")

bulk_tcr_dir <- file.path(data_dir, "bulk_tcrseq")

# ---- Parameters --------------------------------------------------------------
min_templates <- 5     # a clonotype must have >= this many templates in at least one timepoint to be tested
alpha         <- 0.01  # BH-adjusted significance threshold

# `pre_timepoint_label` / `post_timepoint_label` must match values produced
# by `parse_sample_name()` below (adjust both to your own visit naming).
pre_timepoint_label  <- "pre"
post_timepoint_label <- "post"

parse_sample_name <- function(filename) {
  base <- str_remove(basename(filename), "\\.tsv$")
  tibble(
    sample  = base,
    patient = str_extract(base, "^[^_]+"),
    visit   = str_extract(base, "(?<=_).+$")
  ) %>%
    mutate(timepoint = case_when(
      str_detect(tolower(visit), "screen|baseline|pre")  ~ pre_timepoint_label,
      str_detect(tolower(visit), "post|infusion|allo")   ~ post_timepoint_label,
      TRUE ~ NA_character_
    ))
}

tcr_files <- list.files(bulk_tcr_dir, pattern = "\\.tsv$", full.names = TRUE)
all_tcr <- map_dfr(tcr_files, function(f) {
  meta <- parse_sample_name(f)
  read_tsv(f, show_col_types = FALSE) %>%
    filter(!is.na(amino_acid), frame_type == "In") %>%
    mutate(
      clone_id = paste(amino_acid, v_resolved, j_resolved, sep = "_"),
      patient = meta$patient, visit = meta$visit, timepoint = meta$timepoint
    )
}) %>%
  filter(!is.na(timepoint))

# ---- Per-patient paired binomial test ----------------------------------------
run_pre_vs_post_test <- function(patient_id, all_tcr) {
  message("Testing patient: ", patient_id)
  patient_data <- all_tcr %>% filter(patient == patient_id)

  if (!all(c(pre_timepoint_label, post_timepoint_label) %in% patient_data$timepoint)) {
    warning("Skipping ", patient_id, ": missing pre- or post-timepoint sample")
    return(NULL)
  }

  depth_df <- patient_data %>%
    group_by(timepoint) %>%
    summarise(total_templates = sum(templates, na.rm = TRUE), .groups = "drop")

  total_pre  <- depth_df$total_templates[depth_df$timepoint == pre_timepoint_label]
  total_post <- depth_df$total_templates[depth_df$timepoint == post_timepoint_label]
  expected_post_probability <- total_post / (total_pre + total_post)

  paired_tcr <- patient_data %>%
    group_by(patient, clone_id, amino_acid, v_resolved, j_resolved, timepoint) %>%
    summarise(templates = sum(templates, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = timepoint, values_from = templates, values_fill = 0,
                names_prefix = "templates_") %>%
    mutate(
      total_templates_pre_sample  = total_pre,
      total_templates_post_sample = total_post,
      freq_pre  = .data[[paste0("templates_", pre_timepoint_label)]] / total_pre,
      freq_post = .data[[paste0("templates_", post_timepoint_label)]] / total_post,
      total_clone_templates = .data[[paste0("templates_", pre_timepoint_label)]] +
        .data[[paste0("templates_", post_timepoint_label)]],
      testable = .data[[paste0("templates_", pre_timepoint_label)]] >= min_templates |
        .data[[paste0("templates_", post_timepoint_label)]] >= min_templates,
      log2FC_post_vs_pre = log2((freq_post + 1 / (2 * total_post)) / (freq_pre + 1 / (2 * total_pre)))
    )

  results <- paired_tcr %>%
    rowwise() %>%
    mutate(
      p_value = if (testable) {
        binom.test(
          x = .data[[paste0("templates_", post_timepoint_label)]],
          n = total_clone_templates,
          p = expected_post_probability,
          alternative = "two.sided"
        )$p.value
      } else NA_real_
    ) %>%
    ungroup() %>%
    mutate(
      p_adj_BH_within_patient = p.adjust(p_value, method = "BH"),
      direction = case_when(
        is.na(p_adj_BH_within_patient) ~ "Not tested",
        p_adj_BH_within_patient >= alpha ~ "Not significant",
        p_adj_BH_within_patient < alpha & log2FC_post_vs_pre > 0 ~ "Post > Pre",
        p_adj_BH_within_patient < alpha & log2FC_post_vs_pre < 0 ~ "Pre > Post",
        TRUE ~ "Not significant"
      ),
      significant = p_adj_BH_within_patient < alpha
    )
  results
}

all_patients <- sort(unique(all_tcr$patient))
all_results <- map_dfr(all_patients, run_pre_vs_post_test, all_tcr = all_tcr) %>%
  mutate(p_adj_BH_global = p.adjust(p_value, method = "BH"))

significant_results <- all_results %>%
  filter(significant) %>%
  arrange(patient, direction, p_adj_BH_within_patient)

write.csv(all_results, file.path(csv_dir, paste0("clonal_expansion_all_tested_", project_tag, ".csv")), row.names = FALSE)
write.csv(significant_results, file.path(csv_dir, paste0("clonal_expansion_significant_", project_tag, ".csv")), row.names = FALSE)

message("Tested ", n_distinct(all_results$clone_id), " clonotypes across ", length(all_patients),
        " patients; ", nrow(significant_results), " significant (BH q < ", alpha, ") pre/post changes.")
