## ============================================================================
## 05_multimer_calling_and_merge.R
##
## Antigen-specificity ("multimer"/tetramer) calling from the ADT (Antibody
## Capture) assay, followed by merging all per-sample Seurat objects into one
## combined object for downstream integration and clustering.
##
## This study multiplexes each 10x GEM well with a small panel of pMHC
## multimers (see 01_quality_control/feature_references/*.csv), one of which
## targets the HLA-A*02:01-restricted CMV pp65(495-503) epitope NLVPMVATV
## ("NLV"). Multimer identity is called from the highest CLR-normalised ADT
## signal per cell, then confirmed against a sample-specific background
## cutoff. Four antibody-capture products (`multimer_products`, defined in
## 00_config.R) are processed identically -- CYTCMVA001, CYTCMVA002,
## CYTCMVA003 and CYTCMVA006 in this deposit -- each carrying its own
## cutoff/relabelling table below because staining background differs
## slightly between antibody lots/runs.
##
## Input : output/rds/04_seurat_objects_cellcycle_scored_<project_tag>.rds
##         (optional) output/csv/tcr_chain_qc_per_cell_<project_tag>.csv
##           -- produced by 03_tcr_analysis/01_single_cell_vdj/02_tcr_chain_qc_and_collapse.R
##           If not yet available, this script still runs; TCR columns are
##           simply left blank and can be (re)joined later.
## Output: output/plots/<sample>_ADT_ridge_plot.jpeg / _ADT_vln_plot.jpeg
##         output/rds/05_multimer_merged_<project_tag>.rds
## ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(stringr)
library(matrixStats)
library(Matrix)
library(dplyr)
library(tibble)
library(purrr)
library(tidyverse)

source("00_config.R")

seurat_list <- readRDS(file.path(rds_dir, paste0("04_seurat_objects_cellcycle_scored_", project_tag, ".rds")))

#-------------------------------------------------------------------------#
# 1. Per-sample: normalise ADT, call highest-signal multimer per cell
#-------------------------------------------------------------------------#
for (sample_name in names(seurat_list)) {
  seurat <- seurat_list[[sample_name]]
  message("Processing sample: ", sample_name)

  DefaultAssay(seurat) <- "ADT"

  # Exclude hashtag/HTO features -- multimer calling only uses pMHC reagents
  adt_all <- rownames(seurat[["ADT"]])
  adt_features <- adt_all[!grepl("Hashtag|HTO", adt_all, ignore.case = TRUE)]

  seurat <- NormalizeData(seurat, normalization.method = "CLR", features = adt_features)
  seurat <- ScaleData(seurat, features = adt_features)

  p_adt_ridge <- RidgePlot(
    object = seurat, assay = "ADT", layer = "scale.data",
    features = adt_features, y.max = 5, ncol = 5, group.by = "sample"
  ) + ggtitle(paste("ADT staining ridge plot -", sample_name))
  ggsave(filename = file.path(plot_dir, paste0(sample_name, "_ADT_ridge_plot.jpeg")),
         plot = p_adt_ridge, height = 25, width = 25)

  p_vln <- VlnPlot(seurat, features = adt_features, assay = "ADT", layer = "data",
                    group.by = "sample", pt.size = 0.001, ncol = 5)
  ggsave(plot = p_vln, filename = file.path(plot_dir, paste0(sample_name, "_ADT_vln_plot.jpeg")),
         width = 30, height = 30)

  # Assign multimer identity from the highest CLR-normalised ADT signal
  adt_mat <- as.matrix(GetAssayData(seurat, assay = "ADT", layer = "data")[adt_features, ])
  tetramer_call <- apply(adt_mat, 2, function(x) {
    if (all(x == 0)) "negative" else names(x)[which.max(x)]
  })
  seurat <- AddMetaData(seurat, metadata = tetramer_call, col.name = "tetramer")

  adt_data <- GetAssayData(seurat, assay = "ADT", layer = "data")
  tetramer_clr <- map2_chr(colnames(adt_data), seurat$tetramer, function(cell, tetramer_feature) {
    if (tetramer_feature == "negative") NA_character_ else as.character(adt_data[tetramer_feature, cell])
  }) %>% as.numeric()
  seurat$tetramer_CLR <- tetramer_clr

  seurat_list[[sample_name]] <- seurat
}

#-------------------------------------------------------------------------#
# 2. Apply per-sample, per-multimer background cutoffs
#-------------------------------------------------------------------------#
# Illustrative cutoffs for the four antibody-capture products in this panel.
# Replace with values derived from your own negative-control / FMO staining.
# Every product includes the NLV (CMV pp65 NLVPMVATV) reagent.
cutoff_df <- tribble(
  ~sample_name,        ~tetramer,   ~cutoff,
  multimer_products[1], "VTE*TPR",  1.8,
  multimer_products[1], "YSE*CRV",  2.0,
  multimer_products[1], "NLV*QYD",  2.0,
  multimer_products[1], "ELK",      2.1,
  multimer_products[1], "ELR",      2.2,

  multimer_products[2], "VTE*TPR",  1.8,
  multimer_products[2], "YSE*CRV",  3.0,
  multimer_products[2], "NLV*QYD",  2.0,
  multimer_products[2], "ELK",      3.4,
  multimer_products[2], "ELR",      3.4,

  multimer_products[3], "VTE*TPR",  1.8,
  multimer_products[3], "YSE*CRV",  2.0,
  multimer_products[3], "NLV*QYD",  2.0,
  multimer_products[3], "ELK",      3.4,
  multimer_products[3], "ELR",      3.4,

  multimer_products[4], "NLV*QYD",  1.0
)

stopifnot(setequal(unique(cutoff_df$sample_name), multimer_products))

for (sample_name in names(seurat_list)) {
  message("Applying tetramer thresholds for: ", sample_name)
  seurat <- seurat_list[[sample_name]]

  meta <- seurat@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      sample_name  = sample_name,
      tetramer     = as.character(tetramer),
      tetramer_CLR = as.numeric(tetramer_CLR)
    )

  if (!"tetramer_raw" %in% colnames(meta)) {
    meta <- meta %>% mutate(tetramer_raw = tetramer)
  }

  meta_cut <- meta %>%
    left_join(cutoff_df, by = c("sample_name", "tetramer")) %>%
    mutate(
      final_tetramer = case_when(
        tetramer == "negative" ~ "negative",
        !is.na(cutoff) & !is.na(tetramer_CLR) & tetramer_CLR >= cutoff ~ tetramer,
        !is.na(cutoff) & !is.na(tetramer_CLR) & tetramer_CLR < cutoff  ~ "negative",
        TRUE ~ tetramer
      )
    ) %>%
    select(-sample_name, -cutoff) %>%
    column_to_rownames("cell")

  seurat@meta.data <- meta_cut
  seurat_list[[sample_name]] <- seurat
}

# Sanity check: final multimer call distribution per sample
lapply(seurat_list, function(x) table(x$final_tetramer, useNA = "ifany"))

#-------------------------------------------------------------------------#
# 3. Relabel raw feature IDs to short epitope tags per product
#-------------------------------------------------------------------------#
# Which short label a given raw ADT feature maps to can differ slightly
# between products because panels/lots evolved across runs -- e.g. the NLV
# reagent is "NLV*QYD" on three products and simplifies straight to "NLV" on
# the fourth product where it was run alone against background reagents.
rename_map <- tribble(
  ~sample_name,          ~old_label,   ~new_label,

  multimer_products[1],  "VTE*TPR",    "VTE",
  multimer_products[1],  "YSE*CRV",    "YSE",
  multimer_products[1],  "NLV*QYD",    "NLV",

  multimer_products[2],  "VTE*TPR",    "TPR",
  multimer_products[2],  "YSE*CRV",    "CRV",
  multimer_products[2],  "NLV*QYD",    "negative",

  multimer_products[3],  "VTE*TPR",    "VTE",
  multimer_products[3],  "YSE*CRV",    "YSE",
  multimer_products[3],  "NLV*QYD",    "QYD",
  multimer_products[3],  "ELK",        "negative",
  multimer_products[3],  "ELR",        "negative",

  multimer_products[4],  "NLV*QYD",    "NLV",
  multimer_products[4],  "VTE*TPR",    "negative",
  multimer_products[4],  "YSE*CRV",    "negative",
  multimer_products[4],  "ELK",        "negative",
  multimer_products[4],  "ELR",        "negative"
)

for (sample_name in names(seurat_list)) {
  message("Relabelling final_tetramer for: ", sample_name)
  seurat <- seurat_list[[sample_name]]

  meta <- seurat@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      final_tetramer = trimws(as.character(final_tetramer)),
      sample_name = sample_name
    ) %>%
    left_join(rename_map, by = c("sample_name", "final_tetramer" = "old_label")) %>%
    mutate(final_tetramer = ifelse(!is.na(new_label), new_label, final_tetramer)) %>%
    select(-sample_name, -new_label) %>%
    column_to_rownames("cell")

  seurat@meta.data <- meta
  seurat_list[[sample_name]] <- seurat
}

lapply(seurat_list, function(x) table(x$final_tetramer, useNA = "ifany"))

#-------------------------------------------------------------------------#
# 4. Join previously QC'd single-cell TCR chain calls, if available
#-------------------------------------------------------------------------#
# See 03_tcr_analysis/01_single_cell_vdj/02_tcr_chain_qc_and_collapse.R for
# how this file is produced (one row per cell, one TRA/TRB chain each).
tcr_qc_file <- file.path(csv_dir, paste0("tcr_chain_qc_per_cell_", project_tag, ".csv"))

clean_barcode <- function(x) {
  x <- stringr::str_trim(as.character(x))
  bc <- stringr::str_extract(x, "[ACGT]{10,}-[0-9]+$")
  bc <- ifelse(is.na(bc), x, bc)
  stringr::str_remove(bc, "-[0-9]+$")
}

if (file.exists(tcr_qc_file)) {
  tcr_cols <- c("barcode", "TRB_cdr3", "TRB_v_gene", "TRB_j_gene", "TRB_cdr3_nt",
                "TRA_cdr3", "TRA_v_gene", "TRA_j_gene", "TRA_cdr3_nt")

  tcr_raw <- read_csv(tcr_qc_file, show_col_types = FALSE)
  tcr_meta <- tcr_raw %>%
    select(any_of(tcr_cols)) %>%
    mutate(barcode = as.character(barcode), barcode_clean = clean_barcode(barcode)) %>%
    group_by(barcode_clean) %>%
    summarise(across(any_of(setdiff(tcr_cols, "barcode")), ~ {
      vals <- unique(na.omit(as.character(.x))); vals <- vals[vals != ""]
      if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
    }), .groups = "drop")

  for (sample_name in names(seurat_list)) {
    message("Adding TCR chain metadata to: ", sample_name)
    seurat <- seurat_list[[sample_name]]

    meta <- seurat@meta.data %>% rownames_to_column("cell")
    if (!"barcode" %in% colnames(meta)) meta$barcode <- meta$cell

    meta <- meta %>%
      mutate(barcode_meta_clean = clean_barcode(barcode),
             barcode_cell_clean = clean_barcode(cell))

    overlap_meta <- sum(meta$barcode_meta_clean %in% tcr_meta$barcode_clean)
    overlap_cell <- sum(meta$barcode_cell_clean %in% tcr_meta$barcode_clean)
    meta$barcode_join <- if (overlap_cell > overlap_meta) meta$barcode_cell_clean else meta$barcode_meta_clean

    meta <- meta %>% select(-any_of(setdiff(tcr_cols, "barcode")))

    meta_new <- meta %>%
      left_join(tcr_meta, by = c("barcode_join" = "barcode_clean")) %>%
      select(-barcode_meta_clean, -barcode_cell_clean, -barcode_join) %>%
      column_to_rownames("cell")

    seurat@meta.data <- meta_new
    seurat_list[[sample_name]] <- seurat
  }
} else {
  message("No pre-computed TCR chain QC file found at ", tcr_qc_file,
          " -- skipping TCR metadata join. Run the 03_tcr_analysis/01_single_cell_vdj ",
          "scripts first if you want TRA/TRB columns merged in at this stage.")
}

#-------------------------------------------------------------------------#
# 5. Drop the (large) SCT assay, then merge all samples into one object
#-------------------------------------------------------------------------#
for (sample_name in names(seurat_list)) {
  seurat <- seurat_list[[sample_name]]
  if (DefaultAssay(seurat) == "SCT") {
    DefaultAssay(seurat) <- if ("RNA" %in% Assays(seurat)) "RNA" else Assays(seurat)[1]
  }
  if ("SCT" %in% Assays(seurat)) seurat[["SCT"]] <- NULL
  seurat_list[[sample_name]] <- seurat
}

merged_obj <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = project_tag,
  merge.data = TRUE
)

merged_obj@meta.data <- merged_obj@meta.data %>%
  dplyr::select(-any_of("tetramer")) %>%
  dplyr::select(-dplyr::starts_with("pANN")) %>%
  dplyr::select(-dplyr::starts_with("DF"))

saveRDS(merged_obj, file.path(rds_dir, paste0("05_multimer_merged_", project_tag, ".rds")))
message("All samples multimer-called and merged into a single object.")
