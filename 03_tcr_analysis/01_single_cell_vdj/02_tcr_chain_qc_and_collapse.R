## ============================================================================
## 02_tcr_chain_qc_and_collapse.R
##
## Collapses contig-level VDJ data (potentially several TRA/TRB contigs per
## cell) down to a single, best-supported TRA chain and TRB chain per cell:
##  - cells with 2 TRB (or 2 TRA) contigs are resolved by keeping the chain
##    with the higher UMI count
##  - cells where the competing chains have equal (or missing) UMI counts are
##    dropped, since the dominant chain can't be called with confidence
##  - cells with no TRB contig at all are dropped (TRB CDR3 is used as the
##    primary clonotype-defining chain downstream)
##
## Input : output/csv/vdj_contigs_merged_<project_tag>.csv
##         output/rds/04_seurat_objects_cellcycle_scored_<project_tag>.rds (for sample metadata)
## Output: output/csv/tcr_chain_qc_per_cell_<project_tag>.csv
##           (consumed by 01_quality_control/05_multimer_calling_and_merge.R)
## ============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(Seurat)

source("../../00_config.R")

df <- read.csv(file.path(csv_dir, paste0("vdj_contigs_merged_", project_tag, ".csv")))

######## Step 1: one row per cell, chains labelled TRA1/TRA2/TRB1/TRB2 #########

tcr_columns <- c("cdr3", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt")
meta_columns <- setdiff(
  colnames(df),
  c(tcr_columns, "chain", "barcode", "contig_id", "umis", "reads", "length", "productive")
)

df_labeled <- df %>%
  filter(chain %in% c("TRA", "TRB")) %>%
  mutate(across(all_of(tcr_columns), as.character)) %>%
  group_by(barcode, chain) %>%
  mutate(chain_index = row_number()) %>%
  ungroup() %>%
  mutate(chain_label = paste0(chain, chain_index))

df_wide <- df_labeled %>%
  select(barcode, chain_label, all_of(tcr_columns), umis, reads) %>%
  pivot_wider(
    names_from = chain_label,
    values_from = c(all_of(tcr_columns), umis, reads),
    names_glue = "{chain_label}_{.value}"
  )

meta_df <- df %>%
  select(barcode, all_of(meta_columns)) %>%
  distinct(barcode, .keep_all = TRUE)

df_collapsed <- df_wide %>% left_join(meta_df, by = "barcode")

######## Step 2: resolve cells with two TRB chains by UMI count #########

df_collapsed_resolved <- df_collapsed
trb1_cols <- grep("^TRB1_", colnames(df_collapsed_resolved), value = TRUE)
trb2_cols <- grep("^TRB2_", colnames(df_collapsed_resolved), value = TRUE)

multi_trb_targets <- df_collapsed_resolved %>%
  filter(!is.na(TRB1_cdr3) & !is.na(TRB2_cdr3)) %>%
  select(barcode, TRB1_umis, TRB2_umis)

drop_barcodes <- multi_trb_targets %>%
  filter(is.na(TRB1_umis) | is.na(TRB2_umis) | TRB1_umis == TRB2_umis) %>%
  pull(barcode)

best_trb <- multi_trb_targets %>%
  filter(!(barcode %in% drop_barcodes)) %>%
  mutate(TRB_slot = if_else(TRB1_umis > TRB2_umis, "TRB1_id", "TRB2_id")) %>%
  select(barcode, TRB_slot)

df_collapsed_resolved <- df_collapsed_resolved %>% left_join(best_trb, by = "barcode")

for (i in seq_len(nrow(best_trb))) {
  bc <- best_trb$barcode[i]; slot <- best_trb$TRB_slot[i]
  idx <- which(df_collapsed_resolved$barcode == bc)
  if (length(idx) == 0) next
  if (slot == "TRB1_id") {
    df_collapsed_resolved[idx, trb2_cols] <- NA
  } else {
    df_collapsed_resolved[idx, trb1_cols] <- df_collapsed_resolved[idx, trb2_cols]
    df_collapsed_resolved[idx, trb2_cols] <- NA
  }
}

df_collapsed_resolved <- df_collapsed_resolved %>%
  select(-TRB_slot) %>%
  filter(!barcode %in% drop_barcodes)

# Cells with no TRB chain are dropped -- TRB CDR3 defines the clonotype downstream
df_collapsed_resolved <- df_collapsed_resolved %>%
  filter(!is.na(TRB1_cdr3) | !is.na(TRB2_cdr3))

######## Step 3: resolve cells with two TRA chains by UMI count #########

df_collapsed_resolved2 <- df_collapsed_resolved
tra1_cols <- grep("^TRA1_", colnames(df_collapsed_resolved2), value = TRUE)
tra2_cols <- grep("^TRA2_", colnames(df_collapsed_resolved2), value = TRUE)

multi_tra_targets <- df_collapsed_resolved2 %>%
  filter(!is.na(TRA1_cdr3) & !is.na(TRA2_cdr3)) %>%
  select(barcode, TRA1_umis, TRA2_umis)

drop_tra_barcodes <- multi_tra_targets %>%
  filter(is.na(TRA1_umis) | is.na(TRA2_umis) | TRA1_umis == TRA2_umis) %>%
  pull(barcode)

best_tra <- multi_tra_targets %>%
  filter(!(barcode %in% drop_tra_barcodes)) %>%
  mutate(TRA_slot = if_else(TRA1_umis > TRA2_umis, "TRA1_id", "TRA2_id")) %>%
  select(barcode, TRA_slot)

df_collapsed_resolved2 <- df_collapsed_resolved2 %>% left_join(best_tra, by = "barcode")

for (i in seq_len(nrow(best_tra))) {
  bc <- best_tra$barcode[i]; slot <- best_tra$TRA_slot[i]
  idx <- which(df_collapsed_resolved2$barcode == bc)
  if (length(idx) == 0) next
  if (slot == "TRA1_id") {
    df_collapsed_resolved2[idx, tra2_cols] <- NA
  } else {
    df_collapsed_resolved2[idx, tra1_cols] <- df_collapsed_resolved2[idx, tra2_cols]
    df_collapsed_resolved2[idx, tra2_cols] <- NA
  }
}

df_collapsed_resolved2 <- df_collapsed_resolved2 %>%
  select(-any_of("TRA_slot")) %>%
  filter(!barcode %in% drop_tra_barcodes)

######## Step 4: collapse TRA1/TRA2 -> TRA_*, TRB1/TRB2 -> TRB_* #########

df_final <- df_collapsed_resolved2 %>%
  mutate(
    TRB_cdr3    = coalesce(TRB1_cdr3,    TRB2_cdr3),
    TRB_v_gene  = coalesce(TRB1_v_gene,  TRB2_v_gene),
    TRB_d_gene  = coalesce(TRB1_d_gene,  TRB2_d_gene),
    TRB_j_gene  = coalesce(TRB1_j_gene,  TRB2_j_gene),
    TRB_c_gene  = coalesce(TRB1_c_gene,  TRB2_c_gene),
    TRB_cdr3_nt = coalesce(TRB1_cdr3_nt, TRB2_cdr3_nt),
    TRB_umis    = coalesce(TRB1_umis,    TRB2_umis),
    TRB_reads   = coalesce(TRB1_reads,   TRB2_reads),

    TRA_cdr3    = coalesce(TRA1_cdr3,    TRA2_cdr3),
    TRA_v_gene  = coalesce(TRA1_v_gene,  TRA2_v_gene),
    TRA_d_gene  = coalesce(TRA1_d_gene,  TRA2_d_gene),
    TRA_j_gene  = coalesce(TRA1_j_gene,  TRA2_j_gene),
    TRA_c_gene  = coalesce(TRA1_c_gene,  TRA2_c_gene),
    TRA_cdr3_nt = coalesce(TRA1_cdr3_nt, TRA2_cdr3_nt),
    TRA_umis    = coalesce(TRA1_umis,    TRA2_umis),
    TRA_reads   = coalesce(TRA1_reads,   TRA2_reads)
  ) %>%
  select(-starts_with("TRB1_"), -starts_with("TRB2_"),
         -starts_with("TRA1_"), -starts_with("TRA2_")) %>%
  filter(!is.na(TRB_cdr3)) %>%
  select(-any_of("sample"))

######## Step 5: attach sample metadata and write out #########

ref_obj <- readRDS(file.path(rds_dir, paste0("04_seurat_objects_cellcycle_scored_", project_tag, ".rds")))
if (is.list(ref_obj) && !is(ref_obj, "Seurat")) {
  meta_all <- bind_rows(lapply(names(ref_obj), function(nm) {
    m <- ref_obj[[nm]]@meta.data
    m$barcode <- rownames(m)
    m
  }))
} else {
  meta_all <- ref_obj@meta.data
  meta_all$barcode <- rownames(meta_all)
}

meta_to_add <- meta_all %>%
  select(any_of(c("barcode", "orig.ident", "sample"))) %>%
  distinct(barcode, .keep_all = TRUE)

df_final2 <- df_final %>% left_join(meta_to_add, by = "barcode")

df_reordered <- df_final2 %>%
  select(
    barcode, any_of(c("sample", "orig.ident")),
    TRA_v_gene, TRA_d_gene, TRA_j_gene, TRA_c_gene, TRA_cdr3, TRA_cdr3_nt,
    TRB_v_gene, TRB_d_gene, TRB_j_gene, TRB_c_gene, TRB_cdr3, TRB_cdr3_nt,
    everything()
  ) %>%
  filter(!is.na(barcode))

write.csv(df_reordered, file.path(csv_dir, paste0("tcr_chain_qc_per_cell_", project_tag, ".csv")), row.names = FALSE)
message("Wrote one-row-per-cell TCR chain table with ", nrow(df_reordered), " cells.")
