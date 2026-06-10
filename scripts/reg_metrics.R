# Shared regression fold metrics for the TF-rate regression notebooks (13_1, 15_1, 21).
# All in the model's working (log1p) target space, matching the TabPFN CV.

library(tibble)

reg_fold_metrics <- function(pred, truth) {
  tibble(
    rmse     = sqrt(mean((truth - pred)^2)),
    mae      = mean(abs(truth - pred)),
    spearman = suppressWarnings(cor(truth, pred, method = "spearman"))
  )
}
