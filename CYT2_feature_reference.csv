## ============================================================================
## 04_cell_cycle_scoring.R
##
## Scores each cell for S-phase and G2/M-phase signatures (Seurat's built-in
## cell cycle gene lists) so that cell-cycle effects can later be regressed
## out during integration (see 02_clustering/01_integration_and_clustering.R).
##
## Input : output/rds/03_seurat_objects_doublets_removed_<project_tag>.rds
## Output: output/rds/04_seurat_objects_cellcycle_scored_<project_tag>.rds
## ============================================================================

library(Seurat)

source("00_config.R")

seurat_objects <- readRDS(file.path(rds_dir, paste0("03_seurat_objects_doublets_removed_", project_tag, ".rds")))
print("Loaded Seurat objects:")
print(names(seurat_objects))

print("Performing cell cycle scoring and metadata annotation...")

cc.genes <- Seurat::cc.genes
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]

  print(paste("Processing cell cycle scoring for:", sample_name))
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- NormalizeData(seurat_obj, verbose = TRUE)
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

  seurat_objects[[sample_name]] <- seurat_obj
}

saveRDS(seurat_objects, file.path(rds_dir, paste0("04_seurat_objects_cellcycle_scored_", project_tag, ".rds")))
print("Cell cycle scoring completed. Seurat objects saved.")
