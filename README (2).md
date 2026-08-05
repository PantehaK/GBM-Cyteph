# 02_clustering

Integration, unsupervised clustering, cell-type annotation, and downstream
clonal composition/diversity analysis of the merged single-cell object.

Run in order:

| Script | Purpose | Output |
|---|---|---|
| `01_integration_and_clustering.R` | SCTransform (regress cell cycle + %mito), PCA, graph clustering, UMAP, cluster marker genes | `01_clustered_*.rds` |
| `02_celltype_annotation_and_markers.R` | Map clusters to cell-type labels, marker DotPlots, composition-by-sample bar chart | `02_annotated_*.rds` |
| `03_clonal_composition_and_diversity.R` | Per-cell-type clonal expansion, D50 diversity per cell type x product, clone-sharing chord diagrams | plots only |

`03_clonal_composition_and_diversity.R` groups results by the four
antibody-capture products defined in `../00_config.R`
(`multimer_products`), so every plot in this script always breaks results
down across all four products run in this study.

The output of `02_celltype_annotation_and_markers.R`
(`output/rds/02_annotated_*.rds`) is also the expected input for
`../03_tcr_analysis/01_single_cell_vdj/03_export_paired_clones.R`, which
adds cell-type and multimer-specificity annotation to the final exported
clonotype table.
