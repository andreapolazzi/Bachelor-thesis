# Shared fold creation for every modelling analysis (classification CV, regression
# CV and the SHAP out-of-fold runs). One definition, so the same partitions are used
# everywhere and results are directly comparable.
#
#   make_folds(strata, k, n_repeats, seed) -> integer matrix [n x n_repeats]
#       column r = the fold id (1..k) of every sample in repeat r, stratified on
#       `strata` (e.g. alt_status). Reproduces the manual 5x5 construction used in
#       notebooks 02/10/20 when called as make_folds(alt_status).
#
#   SHAP uses a single 5-fold partition = repeat 1 (column 'rep1') of this matrix, so
#   the SHAP out-of-fold split is consistent with the CV folds rather than a separate
#   createFolds() draw.

make_folds <- function(strata, k = 5, n_repeats = 5, seed = 21) {
  strata <- as.factor(strata)
  n <- length(strata)
  set.seed(seed)
  fold_mat <- matrix(0L, n, n_repeats,
                     dimnames = list(NULL, paste0('rep', seq_len(n_repeats))))
  for (r in seq_len(n_repeats)) {
    fid <- integer(n)
    for (lv in levels(strata)) {                 # ALT-low before ALT-high (level order)
      idx <- which(strata == lv)
      fid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    fold_mat[, r] <- fid
  }
  fold_mat
}

# Long list of {train, test, rep, fold} splits from a fold matrix, for CV loops.
folds_to_splits <- function(fold_mat) {
  splits <- list(); id <- 1
  for (r in seq_len(ncol(fold_mat))) {
    for (f in sort(unique(fold_mat[, r]))) {
      splits[[id]] <- list(train = which(fold_mat[, r] != f),
                           test  = which(fold_mat[, r] == f),
                           rep = r, fold = f)
      id <- id + 1
    }
  }
  splits
}
