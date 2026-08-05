# 03_tcr_analysis / 02_bulk_deep_sequencing

Bulk (deep) TCRB repertoire sequencing analysis, independent of the 10x
single-cell data. Written for immunoSEQ-style exports (one `.tsv` per
sample, columns including `amino_acid`, `frame_type`, `templates`,
`v_resolved`, `j_resolved`) loaded either directly (scripts 01-02) or via
[`immunarch`](https://immunarch.com/) (script 03).

| Script | Purpose | Output |
|---|---|---|
| `01_repertoire_diversity_bootstrap.R` | Rarefy every sample to a common template depth and bootstrap clonotype richness, Shannon entropy, and clonality (1 - Pielou evenness), each with a 95% CI | `repertoire_diversity_bootstrap_*.csv` |
| `02_clonal_expansion_significance_testing.R` | Per-patient, per-clonotype exact binomial test for significant expansion/contraction between two timepoints (e.g. pre- vs post-product infusion), BH-corrected | `clonal_expansion_all_tested_*.csv`, `clonal_expansion_significant_*.csv` |
| `03_immunarch_clone_tracking.R` | Track the frequency of the top significantly-expanded clonotypes per patient across every available visit | `top10_expanded_clone_tracking_*.csv`, per-patient kinetics plots |

## Notes

- All three scripts expect raw repertoire files under `data/bulk_tcrseq/`
  and parse `patient`/`visit` from the filename via a small
  `parse_sample_name()` helper defined at the top of each script --
  update that helper to match your own sample naming convention.
- Rarefaction depth (`rarefaction_depth` in script 01) should be set to at
  or below the lowest total template count observed across your samples,
  so diversity metrics are compared on an equal-depth basis.
