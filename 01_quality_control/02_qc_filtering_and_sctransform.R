## ============================================================================
## 02_qc_filtering_and_sctransform.R
##
## Applies standard single-cell QC thresholds (%mito, %ribo, %Hb, min genes/
## UMIs per cell), masks TCR/BCR V(D)J genes prior to normalisation (so that
## clonally expanded receptor transcripts don't dominate variable feature
## selection / PCA), then runs SCTransform + PCA per sample.
##
## Input : output/rds/01_seurat_objects_raw_<project_tag>.rds
## Output: output/rds/02_seurat_objects_sctransformed_<project_tag>.rds
##         output/csv/02_post_qc_summary_<project_tag>.csv
##         output/plots/<sample>_variable_features_<project_tag>.jpg
##         output/plots/<sample>_qc_violin_<project_tag>.jpg
## ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(stringr)
library(parallel)
library(sctransform)
library(dplyr)

source("00_config.R")

########################## QC filtering function ##########################
qc_filtering <- function(seurat_obj, sample_name, qc_thresholds) {
  cat("Applying QC filtering to:", sample_name, "\n")

  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL", col.name = "percent.ribo")
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]", col.name = "percent.hb")

  seurat_obj <- subset(
    seurat_obj,
    subset = percent.ribo < qc_thresholds$percent.ribo &
      percent.hb < qc_thresholds$percent.hb &
      percent.mt < qc_thresholds$percent.mt &
      nFeature_RNA > qc_thresholds$nFeature_RNA_min &
      nCount_RNA > qc_thresholds$nCount_RNA_min
  )
  seurat_obj
}

########################## Receptor-gene masking + SCTransform ##########################
# TCR/BCR V(D)J genes are removed from the RNA counts matrix before
# normalisation/HVG selection/PCA so that clonotype-driven transcript
# abundance doesn't distort clustering; VDJ information itself is added back
# from CellRanger VDJ output later in the pipeline (see 05_ and 03_tcr_analysis/).
immune_receptor_filtering <- function(seurat_obj, sample_name) {
  cat("Applying immune receptor gene masking to:", sample_name, "\n")
  DefaultAssay(seurat_obj) <- "RNA"

  receptor_genes <- grep("^TR[AB]|^IGL|^IGK|^IGHV", rownames(seurat_obj), value = TRUE)
  allowed_igh_genes <- grep("^IGH[MGDHE]$", rownames(seurat_obj), value = TRUE)
  receptor_genes <- c(receptor_genes, allowed_igh_genes)

  non_receptor_genes <- setdiff(rownames(seurat_obj[["RNA"]]), receptor_genes)
  rna_count_filtered <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[non_receptor_genes, ]

  seurat_obj[["RNA"]] <- CreateAssay5Object(counts = rna_count_filtered)

  options(future.globals.maxSize = 2000 * 1024^2)
  seurat_obj <- SCTransform(
    seurat_obj,
    vars.to.regress = "percent.mt",
    variable.features.n = 2000,
    ncells = 3000,
    return.only.var.genes = FALSE
  )

  variable_features <- VariableFeatures(seurat_obj)
  if (length(variable_features) == 0) {
    cat("No variable features found for", sample_name, "\n")
  } else {
    cat("Number of variable features identified for", sample_name, ":", length(variable_features), "\n")
    tryCatch({
      seurat_obj <- RunPCA(seurat_obj, features = variable_features)
      cat("PCA completed successfully for", sample_name, "\n")
    }, error = function(e) {
      cat("PCA failed for", sample_name, ":", e$message, "\n")
    })
  }

  seurat_obj
}

########################## Plotting helpers ##########################
save_variable_feature_plot <- function(seurat_obj, sample_name, out_dir) {
  top10 <- head(VariableFeatures(seurat_obj), 10)
  cat("Top 10 variable features for", sample_name, ":", paste(top10, collapse = ", "), "\n")

  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

  tryCatch({
    ggsave(filename = file.path(out_dir, paste0(sample_name, "_variable_features_", project_tag, ".jpg")),
           plot = plot2, height = 10, width = 10)
  }, error = function(e) {
    cat("Failed to save variable feature plot for", sample_name, ":", e$message, "\n")
  })
}

save_violin_plot <- function(seurat_obj, sample_name, out_dir) {
  p_violin <- VlnPlot(
    seurat_obj,
    features = c("percent.mt", "percent.ribo", "percent.hb", "nCount_RNA", "nFeature_RNA"),
    group.by = "sample",
    ncol = 2
  )
  ggsave(filename = file.path(out_dir, paste0(sample_name, "_qc_violin_", project_tag, ".jpg")),
         plot = p_violin, height = 10, width = 15)
}

########################## Main ##########################
seurat_objects <- readRDS(file.path(rds_dir, paste0("01_seurat_objects_raw_", project_tag, ".rds")))

# QC thresholds -- tune to your assay/tissue
qc_thresholds <- list(
  percent.mt        = 15,
  percent.ribo      = 60,
  percent.hb        = 0.1,
  nFeature_RNA_min  = 200,
  nCount_RNA_min    = 200
)

cat("Starting QC filtering + SCTransform for", length(seurat_objects), "samples...\n")

for (sample_name in names(seurat_objects)) {
  cat("\nProcessing sample:", sample_name, "\n")

  seurat_obj <- seurat_objects[[sample_name]]
  seurat_obj <- qc_filtering(seurat_obj, sample_name, qc_thresholds)
  seurat_obj <- immune_receptor_filtering(seurat_obj, sample_name)

  save_variable_feature_plot(seurat_obj, sample_name, plot_dir)
  save_violin_plot(seurat_obj, sample_name, plot_dir)

  seurat_objects[[sample_name]] <- seurat_obj
}

saveRDS(seurat_objects, file.path(rds_dir, paste0("02_seurat_objects_sctransformed_", project_tag, ".rds")))
cat("All samples QC-filtered, receptor-masked, normalised, and saved.\n")

# Post-filtering QC summary
qc_summary <- list()
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  qc_summary[[sample_name]] <- tibble(
    sample = sample_name,
    nCells = ncol(seurat_obj),
    median_nCount_RNA   = median(seurat_obj$nCount_RNA),
    median_nFeature_RNA = median(seurat_obj$nFeature_RNA),
    median_log10_UMI    = median(log10(seurat_obj$nCount_RNA + 1))
  )
}
write.csv(bind_rows(qc_summary), file.path(csv_dir, paste0("02_post_qc_summary_", project_tag, ".csv")), row.names = FALSE)
