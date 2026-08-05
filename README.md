## ============================================================================
## 03_doublet_removal.R
##
## Runs DoubletFinder on each sample independently: a quick clustering pass
## is used to estimate the homotypic doublet proportion, pK is chosen via a
## parameter sweep, and predicted doublets are removed.
##
## Input : output/rds/02_seurat_objects_sctransformed_<project_tag>.rds
## Output: output/rds/03_seurat_objects_doublets_removed_<project_tag>.rds
##         output/csv/03_post_doublet_removal_summary_<project_tag>.csv
##         output/plots/<sample>_pK_selection_<project_tag>.jpg
## ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(DoubletFinder)
library(reticulate)
library(dplyr)

source("00_config.R")

# Assumed/expected doublet rate for 10x Chromium loading (adjust to your
# target cell recovery per the 10x doublet-rate table).
expected_doublet_rate <- 0.08

seurat_objects <- readRDS(file.path(rds_dir, paste0("02_seurat_objects_sctransformed_", project_tag, ".rds")))
print("Loaded Seurat objects:")
print(names(seurat_objects))

for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]
  print(paste("Processing DoubletFinder for:", sample_name))

  if (ncol(seurat_obj) < 100) {
    print(paste("Skipping", sample_name, "- not enough cells to run DoubletFinder."))
    next
  }

  DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

  options(future.globals.maxSize = 2000 * 1024^2)

  # pK parameter sweep
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  p_pk <- ggplot(bcmvn, aes(x = as.numeric(as.character(pK)), y = BCmetric)) +
    geom_point() + geom_line() +
    labs(title = paste("pK selection for", sample_name), x = "pK", y = "BCmvn metric") +
    theme_minimal()
  ggsave(file.path(plot_dir, paste0(sample_name, "_pK_selection_", project_tag, ".jpg")),
         plot = p_pk, height = 5, width = 5, dpi = 300)

  pK_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  print(paste("Optimal pK for", sample_name, ":", pK_value))

  # Homotypic doublet proportion + expected doublets
  annotations <- seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(expected_doublet_rate * ncol(seurat_obj))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic.prop))

  seurat_obj <- doubletFinder(seurat_obj,
                               PCs = 1:10,
                               pN = 0.25,
                               pK = pK_value,
                               nExp = nExp_poi_adj,
                               sct = TRUE)

  doublet_column <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(doublet_column) == 0) {
    print(paste("Skipping", sample_name, "- DoubletFinder classification column missing."))
    next
  }

  print(paste("Filtering out doublets from:", sample_name))
  seurat_obj_filtered <- subset(seurat_obj, subset = !!as.name(doublet_column) == "Singlet")
  print(paste("Doublet removal complete for", sample_name, "- retained", ncol(seurat_obj_filtered), "cells."))

  seurat_objects[[sample_name]] <- seurat_obj_filtered
}

print("DoubletFinder processing completed for all Seurat objects!")

# Post-doublet-removal QC summary
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
write.csv(bind_rows(qc_summary),
          file.path(csv_dir, paste0("03_post_doublet_removal_summary_", project_tag, ".csv")),
          row.names = FALSE)

saveRDS(seurat_objects, file.path(rds_dir, paste0("03_seurat_objects_doublets_removed_", project_tag, ".rds")))
