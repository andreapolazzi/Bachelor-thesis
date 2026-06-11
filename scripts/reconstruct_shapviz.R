# Reconstruct shapviz objects from the TabPFN SHAP outputs.
#
# Each dir contains: shap_values.csv, X_test_explained.csv, shap_baseline.txt,
# feature_names.txt, shap_space.txt
#
# IMPORTANT space difference:
#   imputation     -> SHAP in PROBABILITY units; baseline ~0.058 = mean P(ALT-high).
#                     Read as 'feature moves P(ALT-high) by +/- x'.
#   recontextual.  -> SHAP in TabPFN's native logit-like space (NOT probability);
#                     baseline ~4.37 is the empty-coalition reference. Only the
#                     ordering / relative magnitudes are interpretable; do NOT
#                     read values or sigmoid(baseline) as a probability.
#

library(shapviz)
library(here)

load_shap <- function(dir) {
  shap <- as.matrix(read.csv(file.path(dir, 'shap_values.csv'), check.names = FALSE))
  X    <- read.csv(file.path(dir, 'X_test_explained.csv'), check.names = FALSE)
  baseline <- as.numeric(readLines(file.path(dir, 'shap_baseline.txt')))
  space    <- readLines(file.path(dir, 'shap_space.txt'))
  stopifnot(ncol(shap) == ncol(X), nrow(shap) == nrow(X))
  message(sprintf('Loaded %s: %d rows x %d features | space=%s | baseline=%.4f',
                  dir, nrow(shap), ncol(shap), space, baseline))
  shapviz(object = shap, X = X, baseline = baseline)
}

# Load the row-aligned metadata (donor_id, cancer_type, cancer_group, alt_status)
# written alongside the OOF outputs. Returns NULL if meta.csv is absent.
load_meta <- function(dir) {
  path <- file.path(dir, 'meta.csv')
  if (!file.exists(path)) {
    message('No meta.csv in ', dir)
    return(NULL)
  }
  read.csv(path, check.names = FALSE)
}

# Rank-average normalization. Per-feature rank scaled to [0,1] (NA-safe), and a
# helper that returns a shapviz object whose X is replaced by these ranks so that
# beeswarm colors and dependence axes are robust to outliers and differing feature
# scales. SHAP values and baseline are untouched. Shared by the ALT-status and
# regression (tTF/bTF) SHAP chunks in notebook 10.
rank01 <- function(x) {
  n <- sum(!is.na(x))
  (rank(x, na.last = 'keep', ties.method = 'average') - 1) / (n - 1)
}

recolor_by_rank <- function(sv) {
  X_rank <- as.data.frame(lapply(get_feature_values(sv), rank01))
  shapviz(object   = get_shap_values(sv),
          X        = X_rank,
          baseline = get_baseline(sv))
}

# Dependence plot with an extra discrete *shape* aesthetic (e.g. ALT status) on
# top of the usual continuous color. `group` must be row-aligned with the
# shapviz object (same OOF order), e.g. meta$alt_status. sv_dependence() exposes
# only a color aesthetic, so we attach `group` to the plot's data and map shape
# onto the existing point layer (which pins no fixed shape).
sv_dependence_shape <- function(sv, v, color_var, group,
                                group_name = 'alt_status',
                                shapes = c('ALT-high' = 17, 'ALT-low' = 16)) {
  p <- sv_dependence(sv, v = v, color_var = color_var)
  p$data[[group_name]] <- factor(group)
  p +
    ggplot2::aes(shape = .data[[group_name]]) +
    ggplot2::scale_shape_manual(values = shapes, name = group_name)
}

# Dependence plot that highlights the top-N most frequent levels of a categorical
# `group` (e.g. cancer_type) and mutes everything else into a single gray 'Other'
# class drawn smaller and more transparent, underneath the highlighted points.
# `group` must be row-aligned with the shapviz object (same OOF order).
# x-axis is the feature value (rank if `sv` was rank-recolored), y is SHAP(v).
sv_dependence_top <- function(sv, v, group, top_n = 6,
                              group_name = 'cancer_type', other_label = 'Other',
                              other_color = 'grey80',
                              size_top = 1.7, size_other = 0.9,
                              alpha_top = 0.9, alpha_other = 0.25,
                              shape_group = NULL, shape_name = 'alt_status',
                              shapes = c('ALT-high' = 17, 'ALT-low' = 16)) {
  X <- get_feature_values(sv)
  S <- get_shap_values(sv)
  stopifnot(length(group) == nrow(S), v %in% colnames(S))

  group <- as.character(group)
  freq  <- sort(table(group), decreasing = TRUE)
  top   <- names(freq)[seq_len(min(top_n, length(freq)))]
  grp   <- factor(ifelse(group %in% top, group, other_label),
                  levels = c(sort(top), other_label))

  df <- data.frame(x = X[[v]], shap = S[, v], grp = grp,
                   is_top = grp != other_label)

  pal <- stats::setNames(scales::hue_pal()(length(top)), sort(top))
  pal[other_label] <- other_color

  # optional discrete shape aesthetic (e.g. alt_status) layered on top of color
  aes_pts <- if (is.null(shape_group)) {
    ggplot2::aes(x, shap, color = grp, size = is_top, alpha = is_top)
  } else {
    stopifnot(length(shape_group) == nrow(S))
    df$shp <- factor(shape_group)
    ggplot2::aes(x, shap, color = grp, shape = shp, size = is_top, alpha = is_top)
  }

  p <- ggplot2::ggplot(df, aes_pts) +
    ggplot2::geom_point(data = df[!df$is_top, ]) +   # gray 'Other' underneath
    ggplot2::geom_point(data = df[df$is_top, ]) +    # highlighted on top
    ggplot2::scale_color_manual(values = pal, name = group_name) +
    ggplot2::scale_size_manual(values = c(`FALSE` = size_other, `TRUE` = size_top),
                               guide = 'none') +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = alpha_other, `TRUE` = alpha_top),
                                guide = 'none') +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = size_top, alpha = alpha_top, shape = 16))) +
    ggplot2::labs(x = v, y = paste0('SHAP(', v, ')'))

  if (!is.null(shape_group)) {
    p <- p +
      ggplot2::scale_shape_manual(values = shapes, name = shape_name) +
      ggplot2::guides(shape = ggplot2::guide_legend(
        override.aes = list(size = size_top, alpha = alpha_top)))
  }
  p
}

# Dependence plot that highlights mesenchymal-origin samples (the minority class
# of interest) against a muted gray / smaller / more transparent non-mesenchymal
# background. `group` is the row-aligned cancer_group column, with values
# 'Mesenchymal_origin' / 'Non_Mesenchymal_origin' (see scripts/process_raw_data.R).
# x-axis is the feature value (rank if `sv` was rank-recolored), y is SHAP(v).
sv_dependence_mes <- function(sv, v, group,
                              highlight = 'Mesenchymal_origin',
                              group_name = 'origin',
                              hi_color = 'firebrick', other_color = 'grey80',
                              size_hi = 1.8, size_other = 0.9,
                              alpha_hi = 0.9, alpha_other = 0.25,
                              hi_label = 'Mesenchymal',
                              other_label = 'Non-mesenchymal',
                              shape_group = NULL, shape_name = 'alt_status',
                              shapes = c('ALT-high' = 17, 'ALT-low' = 16)) {
  X <- get_feature_values(sv)
  S <- get_shap_values(sv)
  stopifnot(length(group) == nrow(S), v %in% colnames(S))

  is_hi <- as.character(group) %in% highlight
  grp   <- factor(ifelse(is_hi, hi_label, other_label),
                  levels = c(other_label, hi_label))
  df    <- data.frame(x = X[[v]], shap = S[, v], grp = grp, is_hi = is_hi)

  pal   <- stats::setNames(c(other_color, hi_color), c(other_label, hi_label))
  sizes <- stats::setNames(c(size_other, size_hi),  c(other_label, hi_label))
  alphs <- stats::setNames(c(alpha_other, alpha_hi), c(other_label, hi_label))

  # optional discrete shape aesthetic (e.g. alt_status) layered on top of color
  aes_pts <- if (is.null(shape_group)) {
    ggplot2::aes(x, shap, color = grp, size = grp, alpha = grp)
  } else {
    stopifnot(length(shape_group) == nrow(S))
    df$shp <- factor(shape_group)
    ggplot2::aes(x, shap, color = grp, shape = shp, size = grp, alpha = grp)
  }

  p <- ggplot2::ggplot(df, aes_pts) +
    ggplot2::geom_point(data = df[!df$is_hi, ]) +   # muted non-mesenchymal underneath
    ggplot2::geom_point(data = df[df$is_hi, ]) +    # mesenchymal highlighted on top
    ggplot2::scale_color_manual(values = pal,   name = group_name) +
    ggplot2::scale_size_manual(values = sizes,  guide = 'none') +
    ggplot2::scale_alpha_manual(values = alphs, guide = 'none') +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = size_hi, alpha = alpha_hi, shape = 16))) +
    ggplot2::labs(x = v, y = paste0('SHAP(', v, ')'))

  if (!is.null(shape_group)) {
    p <- p +
      ggplot2::scale_shape_manual(values = shapes, name = shape_name) +
      ggplot2::guides(shape = ggplot2::guide_legend(
        override.aes = list(size = size_hi, alpha = alpha_hi)))
  }
  p
}

sv_imp <- load_shap(here('outputs', 'shap_output', 'output_imputation_full'))   # probability space
sv_rec <- load_shap(here('outputs', 'shap_output', 'output_recontext_full'))    # logit-like space

# --- 5-fold out-of-fold, blood-inclusive (full2) ------------------------
# Each row is explained out-of-fold; baseline used here is the mean of the
# per-row (per-fold) baselines. Per-row baselines are in baselines.csv if you
# need exact per-sample reconstruction.
sv_imp2 <- load_shap(here('outputs', 'shap_output', 'output_imputation_full2'))  # probability
sv_rec2 <- load_shap(here('outputs', 'shap_output', 'output_recontext_full2'))   # logit-like

# Row-aligned labels for coloring/faceting plots (same order as sv_imp2$X rows):
meta2 <- load_meta(here('outputs', 'shap_output', 'output_imputation_full2'))
# e.g. sv_dependence(sv_imp2, v = 'tf_blood_rate', color_var = meta2$alt_status)

# --- Regression OOF SHAP (tTF = tf_primary_rate, bTF = tf_blood_rate) ----
# Values are in TARGET (log1p) space; baseline = mean of per-row baselines.
# imputation is additive & self-verified; recontext (_rec) is the author-recommended
# cross-check in the model's native target space. Guarded so this file still sources
# before the regression job has produced its outputs.
load_shap_if <- function(dir) if (file.exists(file.path(dir, 'shap_values.csv'))) load_shap(dir) else NULL

sv_ttf_imp <- load_shap_if(here('outputs', 'shap_output_ttf_imp'))
sv_ttf_rec <- load_shap_if(here('outputs', 'shap_output_ttf_rec'))
meta_ttf   <- load_meta(here('outputs', 'shap_output_ttf_imp'))

sv_btf_imp <- load_shap_if(here('outputs', 'shap_output_btf_imp'))
sv_btf_rec <- load_shap_if(here('outputs', 'shap_output_btf_rec'))
meta_btf   <- load_meta(here('outputs', 'shap_output_btf_imp'))
