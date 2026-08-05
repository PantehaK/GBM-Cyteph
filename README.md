# scRNA-seq / TCR-seq analysis code

Analysis code accompanying manuscript 
"Phase I Clinical Trial of Allogeneic Cytomegalovirus-specific T Cells in Combination with Pembrolizumab for Recurrent Glioblastoma".
This repository processes 10x Chromium Single Cell 5' Gene Expression +
Antibody Capture (Hashtag + peptide-MHC multimer) + V(D)J (TCR) data, and
separately, deep (bulk) TCRB repertoire sequencing data, from raw
CellRanger/sequencing output through to the figures and tables reported in
the manuscript.

## Repository structure

```
00_config.R                    Shared paths/settings sourced by every script
01_quality_control/            Per-sample QC, doublet removal, multimer calling, merge
  feature_references/            CellRanger antibody-capture feature-reference CSVs
02_clustering/                 Integration, clustering, cell-type annotation, clonal diversity
03_tcr_analysis/
  01_single_cell_vdj/            10x V(D)J contig processing and clonotype export
  02_bulk_deep_sequencing/       Deep TCRB repertoire diversity/tracking (immunarch)
```

Each folder has its own `README.md` with a script-by-script description.
Scripts are numbered in the order they're meant to be run *within* each
folder; see "Execution order" below for how the folders relate to each
other end-to-end.

## Execution order

```
01_quality_control/01_load_and_create_seurat_objects.R
01_quality_control/02_qc_filtering_and_sctransform.R
01_quality_control/03_doublet_removal.R
01_quality_control/04_cell_cycle_scoring.R
   |
   +--> 03_tcr_analysis/01_single_cell_vdj/01_vdj_merge.R
   +--> 03_tcr_analysis/01_single_cell_vdj/02_tcr_chain_qc_and_collapse.R
   |
01_quality_control/05_multimer_calling_and_merge.R   (joins TCR chain calls in)
   |
02_clustering/01_integration_and_clustering.R
02_clustering/02_celltype_annotation_and_markers.R
02_clustering/03_clonal_composition_and_diversity.R
   |
   +--> 03_tcr_analysis/01_single_cell_vdj/03_export_paired_clones.R

03_tcr_analysis/02_bulk_deep_sequencing/*.R   (independent of the above;
                                                 uses its own deep-sequencing input data)
```

## Multiplexed antibody-capture panel (4 products, incl. CMV-specific barcoded multimers)

Each 10x GEM well was co-stained with a shared panel of sample-demultiplexing
Hashtag antibodies and peptide-MHC multimers, split across two
feature-reference CSVs (`01_quality_control/feature_references/`). Four
antibody-capture products were processed with this pipeline
(`multimer_products` in `00_config.R`): **CYTCMVA001, CYTCMVA002,
CYTCMVA003, CYTCMVA006**. Multimer
specificity is called per cell in
`01_quality_control/05_multimer_calling_and_merge.R` and propagated through
clustering and clonotype export. See `01_quality_control/README.md` for
details of the calling/cutoff logic.

## Environment

Analyses were run in R (>= 4.3) with the following key packages:

- `Seurat` (v5, Assay5/`CreateAssay5Object`), `SeuratObject`
- `sctransform`, `DoubletFinder`, `harmony`
- `immunarch`
- `tidyverse` (`dplyr`, `tidyr`, `readr`, `stringr`, `purrr`, `tibble`, `ggplot2`)
- `patchwork`, `cowplot`, `ggpubr`, `ggrepel`, `ggalluvial`, `circlize`, `scales`
- `matrixStats`, `Matrix`
- `openxlsx`, `readxl`

## Data availability

Raw sequencing data (10x Chromium Gene Expression/V(D)J/Antibody Capture
FASTQs) will be available to request from the European Genome-Phenome 
Archive and accession numbers will be updated after the manuscript is published. 
Processed counts, features, barcodes and final Seurat Object can be requested from
Zenodo under the accession number: 10.5281/zenodo.21800156
repository contains analysis code only; no patient-identifiable data are
included. Paths in `00_config.R` are placeholders (`data/`, `output/`) --
point `data_dir` at your own local copy of the raw/processed data before
running.
