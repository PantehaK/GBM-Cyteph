## ============================================================================
## 01_vdj_merge.R
##
## Scans CellRanger `multi` VDJ-T outputs for every sample present in the
## merged GEX Seurat object, reads each sample's
## `filtered_contig_annotations.csv`, tags contigs with sample/batch, and
## concatenates everything into a single long-format VDJ table (one row per
## contig; cells with paired TRA+TRB have two rows).
##
## Input : output/rds/04_seurat_objects_cellcycle_scored_<project_tag>.rds
##           (only used to get the list of valid sample names / cells; any
##            merged object with `orig.ident` populated will work)
##         data/<GEM_well>/<GEM_well>/outs/per_sample_outs/<sample>/vdj_t/
## Output: output/csv/vdj_contigs_merged_<project_tag>.csv
## ============================================================================

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(rlang)
library(tibble)

source("../../00_config.R")

# Any merged/annotated Seurat object with `orig.ident` populated works here;
# we only need it for the set of valid sample names.
ref_obj <- readRDS(file.path(rds_dir, paste0("04_seurat_objects_cellcycle_scored_", project_tag, ".rds")))
sample_names <- if (is.list(ref_obj) && !is(ref_obj, "Seurat")) names(ref_obj) else unique(ref_obj$orig.ident)

extract_batch_from_path <- function(vdj_dir, base_dir) {
  vdj_dir  <- normalizePath(vdj_dir, winslash = "/", mustWork = FALSE)
  base_dir <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  rel <- sub(paste0("^", base_dir, "/?"), "", vdj_dir)
  batch <- strsplit(rel, "/", fixed = TRUE)[[1]][1]
  if (is.na(batch) || batch == "") batch <- NA_character_
  batch
}

message("Scanning for VDJ-T outputs under: ", data_dir)

top_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
gem_dirs <- top_dirs[grepl(gem_well_prefix, basename(top_dirs))]

vdj_entries <- list()

for (gem_dir in gem_dirs) {
  # Expected layout: <GEM_well>/<GEM_well>/outs/per_sample_outs/<sample>/vdj_t
  vdj_t_dirs <- Sys.glob(file.path(gem_dir, "*", "outs", "per_sample_outs", "*", "vdj_t"))
  if (length(vdj_t_dirs) == 0) {
    message("  No vdj_t folders in: ", gem_dir)
    next
  }

  for (vdj_dir in vdj_t_dirs) {
    sample_name <- basename(dirname(vdj_dir))
    if (!sample_name %in% sample_names) next

    vdj_file <- file.path(vdj_dir, "filtered_contig_annotations.csv")
    if (!file.exists(vdj_file)) {
      message("  Missing VDJ file for sample ", sample_name, " at: ", vdj_file)
      next
    }

    dat <- tryCatch(readr::read_csv(vdj_file, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0) {
      message("  Empty VDJ file for sample ", sample_name, ", skipping.")
      next
    }

    batch <- extract_batch_from_path(vdj_dir, data_dir)

    dat <- dat %>%
      dplyr::mutate(
        barcode     = paste0(sample_name, "_", .data$barcode),
        id          = sample_name,
        batch       = batch,
        source_base = data_dir
      )

    key <- paste(sample_name, vdj_dir, sep = "||")
    vdj_entries[[key]] <- dat
    message("  Added ", nrow(dat), " rows [sample: ", sample_name, ", batch: ", batch, "]")
  }
}

final_merged_vdj <- dplyr::bind_rows(vdj_entries)

write.csv(final_merged_vdj, file.path(csv_dir, paste0("vdj_contigs_merged_", project_tag, ".csv")), row.names = FALSE)
message("Wrote merged VDJ contig table with ", nrow(final_merged_vdj), " rows.")
