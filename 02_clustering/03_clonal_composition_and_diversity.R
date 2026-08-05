## ============================================================================
## 03_clonal_composition_and_diversity.R
##
## Downstream single-cell TCR clonal analyses layered on top of the
## annotated clustering result: clonal composition per cell type, D50 clonal
## diversity per cell type x product, and a clone-sharing chord diagram
## between cell types within each of the four antibody-capture products
## (`multimer_products`; see 00_config.R and 01_quality_control/README.md).
##
## Input : output/rds/02_annotated_<project_tag>.rds
## Output: output/plots/03_clonal_composition_<project_tag>.png
##         output/plots/03_D50_percent_per_celltype_per_product_<project_tag>.png
##         output/plots/03_clone_sharing_chords/<product>_<project_tag>.pdf
## ============================================================================

library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(scales)
library(Seurat)
library(circlize)
library(grid)

source("../00_config.R")

obj <- readRDS(file.path(rds_dir, paste0("02_annotated_", project_tag, ".rds")))

celltype_order <- c(
  "Memory CD8+/CD4+", "IFN-stimulated CD8+", "Proliferating CD8+",
  "Cycling CD8+", "Terminal effector CD8+", "NK"
)
celltype_cols <- c(
  "Memory CD8+/CD4+" = "#2166AC", "IFN-stimulated CD8+" = "#67A9CF",
  "Proliferating CD8+" = "#FDB863", "Cycling CD8+" = "#F46D43",
  "Terminal effector CD8+" = "#D73027", "NK" = "#8E0152"
)

# One colour per antibody-capture product (there are always exactly four in
# this panel design -- see multimer_products in 00_config.R).
product_cols <- setNames(
  c("#2C7FB8", "#FDB863", "#F46D43", "#8E0152")[seq_along(multimer_products)],
  multimer_products
)

## ---- 1. Clonal composition per cell type ------------------------------------
clone_comp <- obj@meta.data %>%
  filter(!is.na(TRB_cdr3), TRB_cdr3 != "", !is.na(celltype)) %>%
  mutate(TRB_clone_key = paste(TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "_")) %>%
  count(celltype, TRB_clone_key, name = "n_cells") %>%
  group_by(celltype) %>%
  mutate(clone_rank = rank(-n_cells, ties.method = "first"),
         is_expanded = n_cells > 1) %>%
  summarise(
    n_cells_total  = sum(n_cells),
    n_clones       = n(),
    n_expanded_clones = sum(is_expanded),
    pct_cells_in_expanded_clones = sum(n_cells[is_expanded]) / sum(n_cells) * 100,
    .groups = "drop"
  )

p_clonal <- ggplot(clone_comp, aes(x = celltype, y = pct_cells_in_expanded_clones, fill = celltype)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_cols) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "% cells in expanded (n>1) TRB clones", title = "Clonal composition by cell type") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(plot = p_clonal, height = 6, width = 8, dpi = 200,
       filename = file.path(plot_dir, paste0("03_clonal_composition_", project_tag, ".png")), bg = "transparent")

## ---- 2. D50 clonal diversity per cell type x product ------------------------
# D50 = the minimum percentage of unique clones that account for 50% of
# cells -- a lower D50 indicates a more clonally skewed (less diverse)
# repertoire; a higher D50 indicates a more polyclonal repertoire.
calc_d50_percent <- function(clones) {
  clones <- clones[!is.na(clones) & clones != ""]
  if (length(clones) == 0) return(NA_real_)
  clone_counts <- sort(table(clones), decreasing = TRUE)
  total_cells  <- sum(clone_counts)
  total_clones <- length(clone_counts)
  clones_needed_for_50 <- which(cumsum(clone_counts) >= 0.5 * total_cells)[1]
  (clones_needed_for_50 / total_clones) * 100
}

meta_d50 <- obj@meta.data %>%
  as.data.frame() %>%
  mutate(
    TRB_clone_key = ifelse(!is.na(TRB_cdr3) & TRB_cdr3 != "",
                            paste(TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "_"), NA)
  ) %>%
  filter(!is.na(celltype), sample %in% multimer_products, !is.na(TRB_clone_key)) %>%
  mutate(celltype = factor(celltype, levels = celltype_order[celltype_order %in% unique(celltype)]),
         sample = factor(sample, levels = multimer_products))

d50_df <- meta_d50 %>%
  group_by(celltype, sample) %>%
  summarise(n_cells = n(), n_clones = n_distinct(TRB_clone_key),
            D50_percent = calc_d50_percent(TRB_clone_key), .groups = "drop")

p_d50 <- ggplot(d50_df, aes(x = celltype, y = D50_percent, fill = sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = product_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0)) +
  labs(x = NULL, y = "D50%", fill = "Product", title = "D50% per cell type and product") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_text(face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(plot = p_d50, height = 4.5, width = 8, dpi = 300,
       filename = file.path(plot_dir, paste0("03_D50_percent_per_celltype_per_product_", project_tag, ".png")),
       bg = "transparent")

## ---- 3. Clone-sharing chord diagrams (per product) ---------------------------
# For each of the four products, visualise how many TRB clonotypes are
# shared between pairs of cell types.
chord_out_dir <- file.path(plot_dir, "03_clone_sharing_chords")
dir.create(chord_out_dir, recursive = TRUE, showWarnings = FALSE)

chord_celltype_order <- setdiff(celltype_order, "NK")  # NK cells lack a clonal TCR repertoire

meta_chord <- obj@meta.data %>%
  as.data.frame() %>%
  mutate(
    TRB_clone_key = ifelse(
      !is.na(TRB_cdr3) & TRB_cdr3 != "" & !is.na(TRB_v_gene) & !is.na(TRB_j_gene),
      paste(TRB_v_gene, TRB_j_gene, TRB_cdr3, sep = "__"), NA_character_
    )
  ) %>%
  filter(!is.na(sample), !is.na(celltype), celltype %in% chord_celltype_order, !is.na(TRB_clone_key)) %>%
  distinct(sample, TRB_clone_key, celltype)

make_chord_matrix <- function(product_id, meta_chord, celltype_order) {
  tcr_df <- meta_chord %>%
    filter(sample == product_id) %>%
    mutate(celltype = factor(celltype, levels = celltype_order)) %>%
    distinct(TRB_clone_key, celltype)

  if (nrow(tcr_df) == 0) return(NULL)

  shared <- tcr_df %>%
    inner_join(tcr_df, by = "TRB_clone_key", relationship = "many-to-many") %>%
    filter(as.integer(celltype.x) < as.integer(celltype.y)) %>%
    count(celltype.x, celltype.y, name = "n_shared_clones")

  if (nrow(shared) == 0) return(NULL)
  shared
}

for (product_id in multimer_products) {
  shared <- make_chord_matrix(product_id, meta_chord, chord_celltype_order)
  if (is.null(shared)) {
    message("No shared TRB clones between cell types for product: ", product_id)
    next
  }

  pdf(file.path(chord_out_dir, paste0(product_id, "_", project_tag, ".pdf")), width = 6, height = 6)
  circos.clear()
  chordDiagram(
    as.data.frame(shared),
    grid.col = celltype_cols[chord_celltype_order],
    annotationTrack = "grid",
    preAllocateTracks = 1
  )
  title(paste("TRB clone sharing between cell types -", product_id))
  circos.clear()
  dev.off()
}

message("Clonal composition and diversity analyses complete.")
