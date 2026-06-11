#!/usr/bin/env Rscript
# Generate the TabPFN-CV input dirs for the classification + tTF + bTF regression
# comparisons, using the shared 5x5 stratified folds (scripts/make_folds.R). Mirrors
# the export chunks of notebooks 20 and 21 exactly, so TabPFN trains on the same folds
# and feature encodings as the in-R models.
#
# Run locally (where R + the data live):
#   Rscript scripts/export_cv_inputs.R
# then sync the three dirs to the chiron $BASE next to tabpfn_cv_oof.py and launch
# scripts/run_tabpfn_cv.sh:
#   rsync -a outputs/shap_input_cmp_cv outputs/cmp_reg_ttf_cv outputs/cmp_reg_btf_cv \
#            chiron1:/ngs/iflores/andrea/

suppressMessages({library(readxl); library(dplyr); library(tidyr); library(here)})
source(here('scripts', 'make_folds.R'))

raw <- read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))

write_inputs <- function(dir, X, y, fold_mat) {
  X <- as.data.frame(X)
  d <- here('outputs', dir)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(X,                        file.path(d, 'X_full.csv'))
  readr::write_csv(data.frame(y = y),        file.path(d, 'y_full.csv'))
  readr::write_csv(as.data.frame(fold_mat),  file.path(d, 'folds_5x5.csv'))
  cat(sprintf("  %-18s  n=%d  p=%d  -> outputs/%s\n", dir, nrow(X), ncol(X), dir))
}

cat("Exporting TabPFN-CV inputs (shared 5x5 folds, seed 21):\n")

## --- classification: 14-feature ALT-status cohort (mirrors notebook 20) ---
features14 <- c(
  'tf_primary_rate', 'tf_blood_rate', 'telomere_insertion_rate',
  'telomere_content_log2', 'TERT_FPKM',
  'ATAGGG_singleton_dist', 'CTAGGG_singleton_dist', 'GTAGGG_singleton_dist',
  'TAAGGG_singleton_dist', 'TCAGGG_singleton_dist', 'TGAGGG_singleton_dist',
  'TTCGGG_singleton_dist', 'TTGGGG_singleton_dist', 'TTTGGG_singleton_dist'
)
clf <- raw %>%
  select(all_of(features14), alt_status) %>%
  filter(complete.cases(across(all_of(features14)))) %>%
  mutate(alt_status = factor(alt_status, levels = c('ALT-low', 'ALT-high')))
write_inputs('shap_input_cmp_cv',
             clf %>% select(all_of(features14)),
             as.integer(clf$alt_status == 'ALT-high'),
             make_folds(clf$alt_status))

## --- regression cohorts (mirror notebook 21 / 13_1 / 15_1) ---
prep <- function(filter_outlier) {
  d <- raw %>% drop_na()
  if (filter_outlier) d <- d %>% filter(tf_primary_rate < 10)
  d %>%
    mutate(
      tf_primary_rate         = log1p(tf_primary_rate),
      tf_blood_rate           = log1p(tf_blood_rate),
      telomere_insertion_rate = log1p(telomere_insertion_rate),
      TERT_FPKM               = log1p(TERT_FPKM),
      alt_status              = factor(alt_status, levels = c('ALT-low', 'ALT-high'))
    ) %>%
    select(-donor_id, -cancer_type, -Specimen.Type.Summary, -cancer_group)
}
export_reg <- function(target, filter_outlier, dir) {
  dat <- prep(filter_outlier)
  X   <- model.matrix(reformulate('.', target), dat)[, -1]   # alt_status -> dummy
  write_inputs(dir, X, dat[[target]], make_folds(dat$alt_status))
}
export_reg('tf_primary_rate', TRUE,  'cmp_reg_ttf_cv')
export_reg('tf_blood_rate',   FALSE, 'cmp_reg_btf_cv')

cat("Done. Sync the three dirs to the chiron $BASE next to tabpfn_cv_oof.py.\n")
