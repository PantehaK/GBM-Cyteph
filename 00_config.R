## ============================================================================
## 00_config.R
##
## Shared project configuration, sourced (or copy-pasted) at the top of each
## pipeline script. Centralising paths here means the only thing a new user
## has to edit to run this pipeline on their own data is this one file.
##
## Directory layout assumed by the pipeline (created automatically if missing):
##
##   data/
##     <GEM_well>/<GEM_well>/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix/
##     <GEM_well>/<GEM_well>/outs/per_sample_outs/<sample>/vdj_t/
##   output/
##     rds/      -- intermediate and final Seurat / R objects
##     csv/      -- tabular summaries, exported clone tables
##     plots/    -- QC and clustering figures
##
## "GEM_well" refers to a single 10x Chromium Next GEM run/batch directory,
## e.g. one 10x lane. Sample names inside `data/` should match the
## `sample_id`/`orig.ident` values used throughout the scripts.
## ============================================================================

# ---- Project-level paths (EDIT THESE) --------------------------------------
project_dir <- "."                                   # repo / project root
data_dir    <- file.path(project_dir, "data")         # raw CellRanger outputs
output_dir  <- file.path(project_dir, "output")
rds_dir     <- file.path(output_dir, "rds")
csv_dir     <- file.path(output_dir, "csv")
plot_dir    <- file.path(output_dir, "plots")

for (d in c(data_dir, output_dir, rds_dir, csv_dir, plot_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Project / product naming -----------------------------------------------
# `project_tag` is used purely as a filename suffix so that outputs from
# different projects/products don't collide. Replace with your own study ID.
project_tag <- "study"

# ---- GEM well / batch prefix -------------------------------------------------
# Regex used to identify per-run subfolders under `data_dir`
# (e.g. folders named CYT1, CYT2, ... for four sequencing/staining runs).
gem_well_prefix <- "^CYT"

# ---- Multiplexed antibody / multimer panel -----------------------------------
# This study multiplexes each 10x GEM well with a panel of Hashtag antibodies
# (sample demultiplexing) and peptide-MHC multimers (antigen specificity).
# Four antibody-capture "products" (CYT_A001-CYT_A004 in the manuscript,
# CYTCMVA001/002/003/006 as run in this deposit) were pooled across two
# feature-reference CSVs (see 01_quality_control/feature_references/).
# One multimer in the panel is specific for the HLA-A*02:01-restricted CMV
# pp65 NLVPMVATV ("NLV") epitope -- see 06_multimer_tetramer_calling.R.
multimer_products <- c("CYTCMVA001", "CYTCMVA002", "CYTCMVA003", "CYTCMVA006")

stopifnot(length(multimer_products) == 4)
