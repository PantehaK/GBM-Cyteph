## ============================================================================
## 01_load_and_create_seurat_objects.R
##
## Reads per-sample, filtered CellRanger multi ("Gene Expression" + "Antibody
## Capture") matrices, builds one Seurat object per sample, attaches an ADT
## assay when hashtag/multimer antibody capture data is present, and computes
## baseline per-cell QC metrics (UMI/gene counts, %mito, %ribo, %hemoglobin).
##
## Input : data/<GEM_well>/<GEM_well>/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix/
## Output: output/rds/01_seurat_objects_raw_<project_tag>.rds
##         output/csv/01_sample_qc_summary_<project_tag>.csv
##         output/plots/01_qc_violin_<project_tag>_batch<N>.png
## ============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(stringr)

source("00_config.R")  # defines data_dir, rds_dir, csv_dir, plot_dir, project_tag, gem_well_prefix

#-------------------------------#
# 1. Discover per-sample matrix folders
#-------------------------------#
batch_folders <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE) %>%
  grep(paste0("^.*/", sub("^\\^", "", gem_well_prefix)), ., value = TRUE)

sample_paths <- unlist(lapply(batch_folders, function(batch) {
  list.dirs(batch, recursive = TRUE, full.names = TRUE) %>%
    grep("sample_filtered_feature_bc_matrix$", ., value = TRUE)
}))

if (length(sample_paths) == 0) {
  stop("No sample_filtered_feature_bc_matrix folders found under ", data_dir,
       ". Check `data_dir` and `gem_well_prefix` in 00_config.R.")
}

#-------------------------------#
# 2. Per-sample processing function
#-------------------------------#
sample_processing_simple <- function(data_dir_sample, sample_name) {
  data <- Read10X(data.dir = data_dir_sample)

  # Regular Seurat object (RNA-only placeholder)
  seurat_obj <- CreateSeuratObject(
    counts  = data$`Gene Expression`,
    project = sample_name
  )
  seurat_obj$sample <- sample_name

  # Replace default RNA assay with an Assay5 object
  rna_assay <- CreateAssay5Object(counts = data$`Gene Expression`, min.cells = 3)
  seurat_obj[["RNA"]] <- rna_assay

  # Add ADT assay (hashtags + multimers) if antibody capture data is present
  if ("Antibody Capture" %in% names(data)) {
    adt_assay <- CreateAssay5Object(counts = data$`Antibody Capture`)
    seurat_obj[["ADT"]] <- adt_assay
  }

  # QC metrics
  seurat_obj$log10_UMI <- log10(seurat_obj$nCount_RNA + 1)

  mito_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  ribo_genes <- grep("^RPS|^RPL", rownames(seurat_obj), value = TRUE)
  hb_genes   <- grep("^HB[AB]", rownames(seurat_obj), value = TRUE)

  seurat_obj$percent.mt <- Matrix::colSums(GetAssayData(seurat_obj, assay = "RNA")[mito_genes, , drop = FALSE]) /
    seurat_obj$nCount_RNA * 100
  seurat_obj$percent.ribo <- Matrix::colSums(GetAssayData(seurat_obj, assay = "RNA")[ribo_genes, , drop = FALSE]) /
    seurat_obj$nCount_RNA * 100
  seurat_obj$percent.hb <- Matrix::colSums(GetAssayData(seurat_obj, assay = "RNA")[hb_genes, , drop = FALSE]) /
    seurat_obj$nCount_RNA * 100

  seurat_obj
}

#-------------------------------#
# 3. Loop through samples
#-------------------------------#
seurat_objects   <- list()
qc_metadata_list <- list()
qc_summary       <- list()

for (path in sample_paths) {
  sample_name <- basename(dirname(dirname(path)))
  message("Processing: ", sample_name)

  obj <- sample_processing_simple(path, sample_name)
  obj$sample_id <- sample_name
  seurat_objects[[sample_name]] <- obj

  qc_metadata_list[[sample_name]] <- obj@meta.data %>%
    dplyr::select(nCount_RNA, nFeature_RNA, log10_UMI, percent.mt, percent.ribo, percent.hb) %>%
    mutate(sample_id = sample_name)

  qc_summary[[sample_name]] <- tibble(
    sample = sample_name,
    nCells = ncol(obj),
    median_nCount_RNA   = median(obj$nCount_RNA),
    median_nFeature_RNA = median(obj$nFeature_RNA),
    median_log10_UMI    = median(log10(obj$nCount_RNA + 1))
  )
}

#-------------------------------#
# 4. Save pre-filtering QC summary
#-------------------------------#
qc_summary_df <- bind_rows(qc_summary)
write.csv(qc_summary_df, file.path(csv_dir, paste0("01_sample_qc_summary_", project_tag, ".csv")), row.names = FALSE)

#-------------------------------#
# 5. Combined QC violin plots (batched, up to 16 samples per figure)
#-------------------------------#
combined_metadata <- bind_rows(qc_metadata_list)

long_df <- pivot_longer(
  combined_metadata,
  cols = c(nCount_RNA, nFeature_RNA, log10_UMI, percent.mt, percent.ribo, percent.hb),
  names_to  = "Metric",
  values_to = "Value"
)

sample_groups <- split(unique(long_df$sample_id), ceiling(seq_along(unique(long_df$sample_id)) / 16))

for (i in seq_along(sample_groups)) {
  group_samples <- sample_groups[[i]]
  df_subset <- long_df %>% filter(sample_id %in% group_samples)

  p <- ggplot(df_subset, aes(x = sample_id, y = Value, fill = sample_id)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "none") +
    labs(title = paste("QC violin plot - samples", min((i - 1) * 16 + 1), "to", min(i * 16)))

  ggsave(
    filename = file.path(plot_dir, paste0("01_qc_violin_", project_tag, "_batch", i, ".png")),
    plot = p, width = 16, height = 10, dpi = 300
  )
}

#-------------------------------#
# 6. Save all Seurat objects
#-------------------------------#
saveRDS(seurat_objects, file.path(rds_dir, paste0("01_seurat_objects_raw_", project_tag, ".rds")))
message("Step 1 complete: raw Seurat objects and QC summary saved for ", length(seurat_objects), " samples.")
