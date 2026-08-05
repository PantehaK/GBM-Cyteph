## ============================================================================
## 03_immunarch_clone_tracking.R
##
## Uses immunarch to load the same bulk repertoire files and tracks the
## trajectory (frequency across all visits) of clonotypes flagged as
## significantly expanded post-product-infusion by
## 02_clonal_expansion_significance_testing.R -- e.g. the top 10
## most-expanded clonotypes per patient -- across every visit that patient
## has, including any visits that weren't part of the pre/post test itself.
##
## Input : output/csv/clonal_expansion_significant_<project_tag>.csv
##         data/bulk_tcrseq/*.tsv (immunarch::repLoad-compatible repertoire files)
## Output: output/csv/top10_expanded_clone_tracking_<project_tag>.csv
##         output/plots/clone_tracking_<patient>_<project_tag>.png
## ============================================================================

library(immunarch)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

source("../../00_config.R")

bulk_tcr_dir <- file.path(data_dir, "bulk_tcrseq")

# ---- 1. Load repertoires with immunarch --------------------------------------
immdata <- repLoad(bulk_tcr_dir)  # expects one file per sample; see immunarch docs for supported formats

repertoire_long <- imap_dfr(immdata$data, function(df, sample_name) {
  df %>%
    mutate(
      sample = sample_name,
      total_templates_visit = sum(Clones, na.rm = TRUE),
      calculated_frequency = Clones / total_templates_visit
    )
})

# Adjust to your own sample naming convention
repertoire_long <- repertoire_long %>%
  mutate(
    patient = str_extract(sample, "^[^_]+"),
    visit   = str_extract(sample, "(?<=_).+$"),
    visit_order = case_when(
      str_detect(tolower(visit), "screen|baseline") ~ 0,
      str_detect(visit, "^V[0-9]+$") ~ as.numeric(str_remove(visit, "^V")),
      TRUE ~ 999
    )
  )

# ---- 2. Visit map per patient -------------------------------------------------
visit_summary_by_patient <- repertoire_long %>%
  distinct(patient, visit, visit_order) %>%
  arrange(patient, visit_order) %>%
  group_by(patient) %>%
  summarise(n_visits = n(), visits_present = paste(visit, collapse = ", "),
            last_visit = visit[which.max(visit_order)], .groups = "drop")
print(visit_summary_by_patient)

# ---- 3. Top 10 significantly-expanded clonotypes per patient -----------------
sig_expanded_file <- file.path(csv_dir, paste0("clonal_expansion_significant_", project_tag, ".csv"))
if (!file.exists(sig_expanded_file)) {
  stop("Run 02_clonal_expansion_significance_testing.R first -- ", sig_expanded_file, " not found.")
}

sig_clones <- read.csv(sig_expanded_file) %>%
  filter(direction == "Post > Pre")

top10_per_patient <- sig_clones %>%
  group_by(patient) %>%
  arrange(desc(log2FC_post_vs_pre), .by_group = TRUE) %>%
  mutate(top10_rank = row_number()) %>%
  slice_head(n = 10) %>%
  ungroup()

# ---- 4. Track those clonotypes' frequency across every visit -----------------
tracking_long <- top10_per_patient %>%
  select(patient, clone_id, CDR3 = amino_acid, top10_rank) %>%
  left_join(
    repertoire_long %>% select(patient, visit, visit_order, CDR3, calculated_frequency),
    by = c("patient", "CDR3"), relationship = "many-to-many"
  ) %>%
  filter(!is.na(visit)) %>%
  arrange(patient, top10_rank, visit_order)

write.csv(tracking_long, file.path(csv_dir, paste0("top10_expanded_clone_tracking_", project_tag, ".csv")), row.names = FALSE)

# ---- 5. Per-patient kinetics plot ---------------------------------------------
for (pid in unique(tracking_long$patient)) {
  df_patient <- tracking_long %>% filter(patient == pid)
  if (nrow(df_patient) == 0) next

  p <- ggplot(df_patient, aes(x = reorder(visit, visit_order), y = calculated_frequency,
                               group = CDR3, colour = factor(top10_rank))) +
    geom_line() + geom_point() +
    labs(x = "Visit", y = "Clonotype frequency", colour = "Rank",
         title = paste("Top 10 post-infusion expanded clonotypes -", pid)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(filename = file.path(plot_dir, paste0("clone_tracking_", pid, "_", project_tag, ".png")),
         plot = p, width = 8, height = 5, dpi = 200)
}

message("Tracked top-10 expanded clonotypes across visits for ", n_distinct(tracking_long$patient), " patients.")
