## ============================================================================
## 02_celltype_annotation_and_markers.R
##
## Maps unsupervised clusters to biologically interpretable cell-type labels
## using canonical marker genes, then visualises marker expression and
## per-sample cell-type composition.
##
## Input : output/rds/01_clustered_<project_tag>.rds
## Output: output/rds/02_annotated_<project_tag>.rds
##         output/plots/02_umap_celltypes_<project_tag>.png
##         output/plots/02_dotplot_celltype_markers_<project_tag>.png
##         output/plots/02_celltype_composition_by_sample_<project_tag>.png
## ============================================================================

library(dplyr)
library(ggplot2)
library(scales)
library(Seurat)

source("../00_config.R")

obj <- readRDS(file.path(rds_dir, paste0("01_clustered_", project_tag, ".rds")))

## ---- Map clusters to cell-type labels --------------------------------------
# Cluster numbers are dataset-specific -- inspect `01_top30_cluster_markers_*.csv`
# and canonical marker DotPlots to assign labels for your own data.
cluster_col <- "seurat_clusters"

cluster_map <- c(
  "0" = "Memory CD8+/CD4+",
  "1" = "Terminal effector CD8+",
  "2" = "Proliferating CD8+",
  "3" = "Cycling CD8+",
  "4" = "NK",
  "5" = "Cycling CD8+",
  "6" = "IFN-stimulated CD8+"
)

obj$celltype <- unname(cluster_map[as.character(obj[[cluster_col, drop = TRUE]])])

celltype_order <- c(
  "Memory CD8+/CD4+",
  "IFN-stimulated CD8+",
  "Proliferating CD8+",
  "Cycling CD8+",
  "Terminal effector CD8+",
  "NK"
)
obj$celltype <- factor(obj$celltype, levels = celltype_order)
obj$celltype_ordered <- obj$celltype

celltype_cols <- c(
  "Memory CD8+/CD4+"       = "#2166AC",
  "IFN-stimulated CD8+"    = "#67A9CF",
  "Proliferating CD8+"     = "#FDB863",
  "Cycling CD8+"           = "#F46D43",
  "Terminal effector CD8+" = "#D73027",
  "NK"                     = "#8E0152"
)

p_umap <- DimPlot(obj, group.by = "celltype", label = FALSE, pt.size = 0.1, alpha = 0.4,
                   label.box = FALSE, stroke.size = 0.4, raster = FALSE, repel = TRUE,
                   cols = celltype_cols, label.size = 4)
ggsave(plot = p_umap, height = 6, width = 8, dpi = 300,
       filename = file.path(plot_dir, paste0("02_umap_celltypes_", project_tag, ".png")), bg = "transparent")

## ---- Canonical marker DotPlot -----------------------------------------------
DefaultAssay(obj) <- "RNA"
genes_to_plot <- c("GZMK", "IL7R", "TYROBP", "CCL4", "IFNG", "GZMB", "IFIT2",
                    "RUNX1", "MKI67", "CD8A", "CD4")
genes_present <- genes_to_plot[genes_to_plot %in% rownames(obj)]

p_dot <- DotPlot(obj, features = genes_present, group.by = "celltype") +
  RotatedAxis() +
  scale_colour_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  scale_y_discrete(limits = celltype_order) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Marker expression by cell type")
ggsave(plot = p_dot, height = 4, width = 8, dpi = 300,
       filename = file.path(plot_dir, paste0("02_dotplot_celltype_markers_", project_tag, ".png")), bg = "transparent")

## ---- Cell-type composition by sample ----------------------------------------
plot_df <- obj@meta.data %>%
  dplyr::select(sample, celltype) %>%
  mutate(celltype = factor(celltype, levels = celltype_order)) %>%
  dplyr::count(sample, celltype, .drop = FALSE) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_stack <- ggplot(plot_df, aes(x = sample, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = celltype_cols, breaks = celltype_order, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.title = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Cell type composition by sample")
ggsave(plot = p_stack, height = 7, width = 8, dpi = 200,
       filename = file.path(plot_dir, paste0("02_celltype_composition_by_sample_", project_tag, ".png")), bg = "transparent")

saveRDS(obj, file.path(rds_dir, paste0("02_annotated_", project_tag, ".rds")))
message("Cell-type annotation complete.")
