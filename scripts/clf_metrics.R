# Shared classification fold metrics for the ALT-status notebooks (10, 20).
#
#   avg_precision(score, y_bin) : area under the precision-recall curve
#                                 (sklearn "average precision"); the honest summary
#                                 metric under heavy class imbalance.
#   fold_metrics(prob, truth)   : one row of per-fold metrics. `truth` is a factor
#                                 with levels c("ALT-low", "ALT-high"); thresholded
#                                 metrics use 0.5 with positive = ALT-high.

library(pROC)
library(tibble)

avg_precision <- function(score, y_bin) {
  o  <- order(score, decreasing = TRUE)
  yb <- y_bin[o]
  tp <- cumsum(yb)
  fp <- cumsum(1 - yb)
  P  <- tp / (tp + fp)
  R  <- tp / sum(yb)
  Rprev <- c(0, head(R, -1))
  sum((R - Rprev) * P)
}

fold_metrics <- function(prob, truth) {
  yb   <- as.integer(truth == "ALT-high")
  ro   <- roc(truth, prob, levels = c("ALT-low", "ALT-high"),
              direction = "<", quiet = TRUE)
  pred <- prob >= 0.5
  tp <- sum(pred  & yb == 1); fp <- sum(pred  & yb == 0)
  tn <- sum(!pred & yb == 0); fn <- sum(!pred & yb == 1)
  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  f1   <- if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_
  tibble(
    auc = as.numeric(pROC::auc(ro)), pr_auc = avg_precision(prob, yb),
    sensitivity = sens, specificity = spec, precision = prec, f1 = f1,
    balanced_accuracy = 0.5 * (sens + spec), n_high = tp + fn, n_low = tn + fp
  )
}
