# 03_tcr_analysis / 01_single_cell_vdj

Processing of 10x Chromium Single Cell V(D)J (TCR) data paired with the
5' Gene Expression + Antibody Capture libraries processed in
`../../01_quality_control/`.

| Script | Purpose | Output |
|---|---|---|
| `01_vdj_merge.R` | Concatenate `filtered_contig_annotations.csv` across all samples/GEM wells | `vdj_contigs_merged_*.csv` |
| `02_tcr_chain_qc_and_collapse.R` | Collapse multi-contig cells to one TRA + one TRB chain per cell (UMI-based tie-breaking) | `tcr_chain_qc_per_cell_*.csv` |
| `03_export_paired_clones.R` | Compute per-sample clonotype frequencies and export a paired TRA/TRB clone table annotated with cell type and multimer specificity | `paired_tcr_clones_*.csv` |

## Execution order relative to the other folders

1. Run `01_vdj_merge.R` and `02_tcr_chain_qc_and_collapse.R` after
   `../../01_quality_control/04_cell_cycle_scoring.R` -- the resulting
   `tcr_chain_qc_per_cell_*.csv` is consumed by
   `../../01_quality_control/05_multimer_calling_and_merge.R` to attach
   TRA/TRB calls to the merged object.
2. Run `03_export_paired_clones.R` last, after clustering
   (`../../02_clustering/02_celltype_annotation_and_markers.R`), since it
   annotates each clonotype with its cell type(s) and multimer specificity.
