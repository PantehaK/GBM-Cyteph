## ============================================================================
## 01_integration_and_clustering.R
##
## Re-normalises the merged, multimer-annotated object (SCTransform,
## regressing out cell-cycle and %mito effects), runs PCA, graph-based
## clustering, and UMAP.
##
## Input : output/rds/05_multimer_merged_<project_tag>.rds
##           (see 01_quality_control/05_multimer_calling_and_merge.R)
## Output: output/rds/01_clustered_<project_tag>.rds
##         output/csv/01_top30_cluster_markers_<project_tag>.csv
##         output/plots/01_umap_clusters_<project_tag>.png
## ============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Seurat)
library(ggplot2)
library(harmony)
library(future)

source("../00_config.R")

obj <- readRDS(file.path(rds_dir, paste0("05_multimer_merged_", project_tag, ".rds")))

# Re-join per-sample layers, then re-split by sample so SCTransform models
# technical variation within each sample before re-merging for integration.
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

obj <- SCTransform(
  obj,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
  return.only.var.genes = TRUE,
  verbose = TRUE
)
DefaultAssay(obj) <- "SCT"

obj <- RunPCA(obj, assay = "SCT")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

## ---- Graph-based clustering + UMAP -----------------------------------------
# Swap RunPCA/FindNeighbors for a Harmony-integrated reduction here (e.g.
# `obj <- RunHarmony(obj, group.by.vars = "sample")`, then
# `reduction = "harmony"` below) if batch effects between GEM wells/products
# are substantial in your dataset.
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = "pca_snn")
obj <- FindClusters(obj, resolution = 0.25, graph.name = "pca_snn", cluster.name = "seurat_clusters")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap")

cluster_cols <- c("#6CC3F4", "#c90076", "#F5A56B", "#8CCE8A", "#B560DD", "#f467ba",
                   "#3A68AE", "#8e7cc3", "#E65757", "#C96BF4", "#DD9560", "#d9a813",
                   "#CC4E4E", "#569CC3")

p_clusters <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.1,
                       alpha = 0.4, label.box = TRUE, stroke.size = 0.4, raster = FALSE,
                       repel = TRUE, cols = cluster_cols, label.size = 4)
ggsave(plot = p_clusters, height = 6, width = 8, dpi = 300,
       filename = file.path(plot_dir, paste0("01_umap_clusters_", project_tag, ".png")),
       bg = "transparent")

## ---- Marker genes per cluster ----------------------------------------------
DefaultAssay(obj) <- "SCT"
prepSCT <- PrepSCTFindMarkers(obj, assay = "SCT")
cluster.markers <- FindAllMarkers(prepSCT, only.pos = TRUE, assay = "SCT", slot = "data")

cluster.markers <- cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5)

# Drop ribosomal/mitochondrial genes -- rarely biologically informative and
# tend to dominate marker lists due to their high, variable expression.
housekeeping_genes <- grep("^RPS|^RPL|^MT-", cluster.markers$gene, value = TRUE)
markers_filtered <- cluster.markers[!(cluster.markers$gene %in% housekeeping_genes), ]

top30 <- markers_filtered %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::ungroup()

write.csv(top30, file.path(csv_dir, paste0("01_top30_cluster_markers_", project_tag, ".csv")))

saveRDS(obj, file.path(rds_dir, paste0("01_clustered_", project_tag, ".rds")))
message("Clustering complete: ", length(unique(obj$seurat_clusters)), " clusters identified.")
