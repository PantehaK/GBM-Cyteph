# 01_quality_control

Per-sample preprocessing of 10x Chromium Single Cell 5' Gene Expression +
Antibody Capture (Hashtag + pMHC multimer) data, from raw filtered matrices
through to a single merged, QC'd, multimer-annotated Seurat object.

Run in order:

| Script | Purpose | Output |
|---|---|---|
| `01_load_and_create_seurat_objects.R` | Read CellRanger `multi` outputs, build one Seurat object per sample, compute %mito/%ribo/%Hb | `01_seurat_objects_raw_*.rds` |
| `02_qc_filtering_and_sctransform.R` | Apply QC thresholds, mask TCR/BCR V(D)J genes, run SCTransform + PCA | `02_seurat_objects_sctransformed_*.rds` |
| `03_doublet_removal.R` | Per-sample DoubletFinder | `03_seurat_objects_doublets_removed_*.rds` |
| `04_cell_cycle_scoring.R` | Score S/G2M phase | `04_seurat_objects_cellcycle_scored_*.rds` |
| `05_multimer_calling_and_merge.R` | Call pMHC multimer specificity per cell, merge all samples | `05_multimer_merged_*.rds` |

## Multiplexed antibody-capture panel

Each 10x GEM well in this study was stained with a shared antibody-capture
panel combining sample-demultiplexing hashtags and peptide-MHC multimers
(Immudex Dextramer-type reagents), split across two feature-reference CSVs
because more antibody-capture features were used than fit one 10x kit
config (`feature_references/CYT1_feature_reference.csv` and
`CYT2_feature_reference.csv`).

Four antibody-capture products were run (`multimer_products` in
`../00_config.R`: `CYTCMVA001`, `CYTCMVA002`, `CYTCMVA003`, `CYTCMVA006`),
each pooling the hashtag(s) for that GEM well with the shared multimer
panel. One multimer in the panel is specific for the HLA-A\*02:01-restricted
CMV pp65(495-503) epitope **NLVPMVATV ("NLV")** -- a widely used
immunodominant CMV epitope. Multimer specificity is called per cell from the
highest CLR-normalised antibody-capture signal, then confirmed against a
sample-specific background cutoff (`05_multimer_calling_and_merge.R`).

Feature reference CSVs follow the standard CellRanger `feature-ref.csv`
format (`id,name,read,pattern,sequence,feature_type`) and are provided here
as-run examples for reproducing the multiplexing scheme.
