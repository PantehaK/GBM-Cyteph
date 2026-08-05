## ============================================================================
## 03_export_paired_clones.R
##
## Defines a clonotype as the combination of TRB CDR3 + V + J gene (paired
## with TRA CDR3 + V + J gene when available), computes each clonotype's
## frequency within its sample, and exports one row per clonotype x sample
## with cell-type composition and multimer/antigen-specificity annotation.
##
## Input : output/rds/<clustering output, e.g. 02_clustering/02_annotated_<project_tag>.rds>
## Output: output/csv/paired_tcr_clones_<project_tag>.csv
## ============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Seurat)

source("../../00_config.R")

# Point this at the final annotated/clustered object
# (see 02_clustering/02_celltype_annotation_and_markers.R)
obj2 <- readRDS(file.path(rds_dir, paste0("02_annotated_", project_tag, ".rds")))

collapse_unique <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

collapse_with_counts <- function(x, sep = "; ") {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  paste0(names(tab), " (", as.integer(tab), " cells)", collapse = sep)
}

dominant_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

# Cells double-called across mutually exclusive multimers (rare, background)
# are flagged as "Mixed" rather than arbitrarily assigned to one specificity.
specificity_value <- function(x) {
  x <- x[!is.na(x) & x != "" & !tolower(x) %in% c("negative", "none")]
  x <- unique(x)
  if (length(x) == 0) return("Negative/none")
  if (length(x) == 1) return(x)
  paste("Mixed:", paste(x, collapse = "; "))
}

meta <- obj2@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    TRB_clone_id = if_else(
      !is.na(TRB_cdr3) & TRB_cdr3 != "" & !is.na(TRB_v_gene) & TRB_v_gene != "" &
        !is.na(TRB_j_gene) & TRB_j_gene != "",
      paste(TRB_cdr3, TRB_v_gene, TRB_j_gene, sep = "_"),
      NA_character_
    ),
    TRA_clone_id = if_else(
      !is.na(TRA_cdr3) & TRA_cdr3 != "" & !is.na(TRA_v_gene) & TRA_v_gene != "" &
        !is.na(TRA_j_gene) & TRA_j_gene != "",
      paste(TRA_cdr3, TRA_v_gene, TRA_j_gene, sep = "_"),
      "No_TRA"
    ),
    paired_clone_id = if_else(!is.na(TRB_clone_id), paste(TRB_clone_id, TRA_clone_id, sep = "__"), NA_character_)
  )

paired_clone_export <- meta %>%
  filter(!is.na(paired_clone_id)) %>%
  group_by(sample) %>%
  mutate(total_paired_TCR_cells_sample = n()) %>%
  group_by(sample, paired_clone_id, TRB_clone_id, TRA_clone_id,
           TRB_cdr3, TRB_v_gene, TRB_j_gene, TRA_cdr3, TRA_v_gene, TRA_j_gene) %>%
  summarise(
    n_cells_paired_clone = n(),
    total_paired_TCR_cells_sample = first(total_paired_TCR_cells_sample),
    paired_clone_prop_sample = n_cells_paired_clone / total_paired_TCR_cells_sample,

    TRB_cdr3_nt = collapse_unique(TRB_cdr3_nt),
    TRA_cdr3_nt = collapse_unique(TRA_cdr3_nt),

    # `final_tetramer` is the multimer-specificity call from
    # 01_quality_control/05_multimer_calling_and_merge.R (includes the NLV
    # CMV-specificity call where present).
    final_tetramer_specificity = specificity_value(final_tetramer),
    final_tetramer_dominant    = dominant_value(final_tetramer),
    final_tetramer_all_calls   = collapse_with_counts(final_tetramer),

    celltypes = collapse_unique(celltype),

    .groups = "drop"
  ) %>%
  arrange(sample, desc(n_cells_paired_clone), TRB_cdr3, TRA_cdr3)

write.csv(paired_clone_export, file.path(csv_dir, paste0("paired_tcr_clones_", project_tag, ".csv")), row.names = FALSE)
message("Wrote paired clonotype table with ", nrow(paired_clone_export), " sample x clonotype rows.")
