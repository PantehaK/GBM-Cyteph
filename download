# 03_tcr_analysis

TCR analyses are split into two independent subfolders because they draw
on two different data types:

- **`01_single_cell_vdj/`** -- paired TRA/TRB clonotype calls from 10x
  Chromium Single Cell V(D)J sequencing, linked back to the transcriptomic
  clusters/cell types and multimer-specificity calls from `../02_clustering/`.
- **`02_bulk_deep_sequencing/`** -- deep (bulk) TCRB repertoire sequencing
  (e.g. Adaptive immunoSEQ), used for high-sensitivity detection and
  statistical tracking of clonotype expansion across longitudinal patient
  samples, analysed with `immunarch` and custom bootstrap/rarefaction code.

The two are complementary: the single-cell data give per-clonotype
phenotype/antigen-specificity (including CMV-NLV specificity), while the
bulk deep-sequencing data give the statistical power to detect rare,
significantly expanded clonotypes across timepoints that single-cell
sampling alone would likely miss.

See each subfolder's own README for script-by-script details.
