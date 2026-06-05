# Reconstruct shapviz objects from the TabPFN SHAP outputs.
#
# Each dir contains: shap_values.csv, X_test_explained.csv, shap_baseline.txt,
# feature_names.txt, shap_space.txt
#
# IMPORTANT space difference:
#   imputation     -> SHAP in PROBABILITY units; baseline ~0.058 = mean P(ALT-high).
#                     Read as "feature moves P(ALT-high) by +/- x".
#   recontextual.  -> SHAP in TabPFN's native logit-like space (NOT probability);
#                     baseline ~4.37 is the empty-coalition reference. Only the
#                     ordering / relative magnitudes are interpretable; do NOT
#                     read values or sigmoid(baseline) as a probability.
#

library(shapviz)
library(here)

load_shap <- function(dir) {
  shap <- as.matrix(read.csv(file.path(dir, "shap_values.csv"), check.names = FALSE))
  X    <- read.csv(file.path(dir, "X_test_explained.csv"), check.names = FALSE)
  baseline <- as.numeric(readLines(file.path(dir, "shap_baseline.txt")))
  space    <- readLines(file.path(dir, "shap_space.txt"))
  stopifnot(ncol(shap) == ncol(X), nrow(shap) == nrow(X))
  message(sprintf("Loaded %s: %d rows x %d features | space=%s | baseline=%.4f",
                  dir, nrow(shap), ncol(shap), space, baseline))
  shapviz(object = shap, X = X, baseline = baseline)
}

# Load the row-aligned metadata (donor_id, cancer_type, cancer_group, alt_status)
# written alongside the OOF outputs. Returns NULL if meta.csv is absent.
load_meta <- function(dir) {
  path <- file.path(dir, "meta.csv")
  if (!file.exists(path)) {
    message("No meta.csv in ", dir)
    return(NULL)
  }
  read.csv(path, check.names = FALSE)
}

sv_imp <- load_shap(here("outputs", "shap_output", "output_imputation_full"))   # probability space
sv_rec <- load_shap(here("outputs", "shap_output", "output_recontext_full"))    # logit-like space

# --- 5-fold out-of-fold, blood-inclusive (full2) ------------------------
# Each row is explained out-of-fold; baseline used here is the mean of the
# per-row (per-fold) baselines. Per-row baselines are in baselines.csv if you
# need exact per-sample reconstruction.
sv_imp2 <- load_shap(here("outputs", "shap_output", "output_imputation_full2"))  # probability
sv_rec2 <- load_shap(here("outputs", "shap_output", "output_recontext_full2"))   # logit-like

# Row-aligned labels for coloring/faceting plots (same order as sv_imp2$X rows):
meta2 <- load_meta(here("outputs", "shap_output", "output_imputation_full2"))
# e.g. sv_dependence(sv_imp2, v = "tf_blood_rate", color_var = meta2$alt_status)
