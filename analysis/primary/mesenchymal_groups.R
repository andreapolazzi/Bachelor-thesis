# =============================================================
# Group comparison: Mesenchymal_origin vs Non-Mesenchymal_origin
# Exploratory only — n_mes = 42, n_nonmes = 774
# =============================================================

library(tidyverse)
library(here)
library(patchwork)
library(effsize)   # for Cliff's delta

data <- readxl::read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))
clean_data <- data %>% 
  drop_na() %>% 
  filter(tf_primary_rate<10)

cancer_type_vector <- clean_data$cancer_type
cancer_group_vector <- clean_data$cancer_group

clean_data_orig <- clean_data %>% 
  mutate(
    tf_blood_rate = log1p(tf_blood_rate),
    telomere_insertion_rate = log1p(telomere_insertion_rate),
    TERT_FPKM = log1p(TERT_FPKM),
    alt_status = factor(
      alt_status,
      levels = c("ALT-low", "ALT-high")
    )
  ) %>%  
  select(-donor_id, -cancer_type, -Specimen.Type.Summary, -cancer_group)

clean_data <- clean_data %>% 
  mutate(
    tf_primary_rate = log1p(tf_primary_rate),
    tf_blood_rate = log1p(tf_blood_rate),
    telomere_insertion_rate = log1p(telomere_insertion_rate),
    TERT_FPKM = log1p(TERT_FPKM),
    alt_status = factor(
      alt_status,
      levels = c("ALT-low", "ALT-high")
    )
  ) %>%  
  select(-donor_id, -cancer_type, -Specimen.Type.Summary, -cancer_group)

# ---- 0. Build the working frame -----------------------------

df_group <- clean_data %>%
  mutate(cancer_group = factor(cancer_group_vector))

# Extract group names ONCE from the data, don't hardcode
grp_levels <- levels(df_group$cancer_group)
stopifnot(length(grp_levels) == 2)  # this analysis assumes 2 groups

# Identify which is which by name pattern (robust to '-' vs '_')
grp_mes <- grp_levels[grepl("(?i)^mesenchymal", grp_levels)]
grp_non <- setdiff(grp_levels, grp_mes)
stopifnot(length(grp_mes) == 1, length(grp_non) == 1)

# Reorder so non-mesenchymal is the reference
df_group$cancer_group <- factor(df_group$cancer_group,
                                levels = c(grp_non, grp_mes))

cat("Groups detected:\n")
cat("  Reference:   ", grp_non, "\n")
cat("  Comparison:  ", grp_mes, "\n\n")
print(table(df_group$cancer_group))

# Colour palette built once, keyed on the actual level names
grp_colors      <- setNames(c("grey60", "firebrick"),    c(grp_non, grp_mes))
grp_colors_dark <- setNames(c("grey40", "firebrick4"),   c(grp_non, grp_mes))
grp_sizes       <- setNames(c(1.2, 2.2),                  c(grp_non, grp_mes))
grp_alphas      <- setNames(c(0.4, 0.9),                  c(grp_non, grp_mes))


# ---- 1. Quantify differences for every numeric feature ------

numeric_vars <- df_group %>%
  select(-cancer_group, -alt_status) %>%
  names()

safe_cliff <- function(x, y) {
  if (length(unique(c(x, y))) < 2) {
    return(list(estimate = 0, magnitude = "negligible",
                note = "all values identical"))
  }
  res <- tryCatch(effsize::cliff.delta(x, y), error = function(e) NULL)
  if (is.null(res) || is.nan(res$estimate)) {
    n_x <- length(x); n_y <- length(y)
    greater <- sum(outer(x, y, ">"))
    less    <- sum(outer(x, y, "<"))
    d <- (greater - less) / (n_x * n_y)
    mag <- cut(abs(d),
               breaks = c(-Inf, 0.147, 0.33, 0.474, Inf),
               labels = c("negligible", "small", "medium", "large"),
               right = FALSE)
    return(list(estimate = d, magnitude = as.character(mag),
                note = "manual fallback"))
  }
  list(estimate = as.numeric(res$estimate),
       magnitude = as.character(res$magnitude),
       note = "")
}

safe_smd <- function(x, y) {
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled_sd
}

group_diff_numeric <- map_dfr(numeric_vars, function(v) {
  x_mes <- df_group[[v]][df_group$cancer_group == grp_mes]
  x_non <- df_group[[v]][df_group$cancer_group == grp_non]
  
  cd <- safe_cliff(x_mes, x_non)
  
  tibble(
    feature       = v,
    n_mes         = length(x_mes),
    n_non         = length(x_non),
    n_unique_mes  = length(unique(x_mes)),
    n_unique_non  = length(unique(x_non)),
    median_mes    = median(x_mes),
    median_non    = median(x_non),
    cliff_delta   = cd$estimate,
    cliff_mag     = cd$magnitude,
    smd           = safe_smd(x_mes, x_non),
    note          = cd$note
  )
}) %>%
  arrange(desc(abs(cliff_delta)))

print(group_diff_numeric, n = Inf)

# ALT-stratified effect sizes (with bootstrap CI) — only in ALT-high subset
df_alt_high <- df_group %>% filter(alt_status == "ALT-high")

cliff_with_ci <- function(x, y, n_boot = 4000) {
  delta_obs <- safe_cliff(x, y)$estimate
  boots <- replicate(n_boot, {
    xb <- sample(x, length(x), replace = TRUE)
    yb <- sample(y, length(y), replace = TRUE)
    safe_cliff(xb, yb)$estimate
  })
  boots <- boots[is.finite(boots)]
  list(
    estimate = delta_obs,
    ci_low   = unname(quantile(boots, 0.025)),
    ci_high  = unname(quantile(boots, 0.975))
  )
}

set.seed(21)
group_diff_alt_high <- map_dfr(numeric_vars, function(v) {
  x_mes <- df_alt_high[[v]][df_alt_high$cancer_group == grp_mes]
  x_non <- df_alt_high[[v]][df_alt_high$cancer_group == grp_non]
  cd_ci <- cliff_with_ci(x_mes, x_non)
  mag <- cut(abs(cd_ci$estimate),
             breaks = c(-Inf, 0.147, 0.33, 0.474, Inf),
             labels = c("negligible", "small", "medium", "large"),
             right = FALSE)
  tibble(
    feature     = v,
    n_mes       = length(x_mes),
    n_non       = length(x_non),
    cliff_delta = cd_ci$estimate,
    ci_low      = cd_ci$ci_low,
    ci_high     = cd_ci$ci_high,
    cliff_mag   = as.character(mag),
    excludes_0  = (cd_ci$ci_low > 0) | (cd_ci$ci_high < 0)
  )
}) %>% arrange(desc(abs(cliff_delta)))

print(group_diff_alt_high, n = Inf)
# ---- 2. Effect-size ranking plot ----------------------------

group_diff_plot <- group_diff_numeric %>%
  filter(!is.na(cliff_delta)) %>%
  mutate(
    cliff_mag = factor(cliff_mag,
                       levels = c("negligible", "small", "medium", "large")),
    feature   = fct_reorder(feature, abs(cliff_delta))
  )

plt_effect_ranking <- group_diff_plot %>%
  ggplot(aes(x = cliff_delta, y = feature, fill = cliff_mag)) +
  geom_col(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.147, 0.147),
             linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = c(-0.33, 0.33),
             linetype = "dotted", color = "grey60") +
  scale_fill_manual(
    values = c(negligible = "grey80", small = "#fdcc8a",
               medium = "#fc8d59", large = "#d7301f"),
    drop = FALSE
  ) +
  labs(
    title    = "Cliff's delta: Mesenchymal vs Non-Mesenchymal",
    subtitle = sprintf("Positive = higher in %s (n=%d vs %d). Dotted lines = negligible/small thresholds.",
                       grp_mes,
                       sum(df_group$cancer_group == grp_mes),
                       sum(df_group$cancer_group == grp_non)),
    x = "Cliff's delta", y = NULL, fill = "Effect size"
  ) +
  theme_bw()

plt_effect_ranking

group_diff_plot <- group_diff_alt_high %>%
  filter(!is.na(cliff_delta)) %>%
  mutate(
    cliff_mag = factor(cliff_mag,
                       levels = c("negligible", "small", "medium", "large")),
    feature   = fct_reorder(feature, abs(cliff_delta))
  )

plt_effect_ranking <- group_diff_plot %>%
  ggplot(aes(x = cliff_delta, y = feature, fill = cliff_mag)) +
  geom_col(alpha = 0.85) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.3, color = "grey25") +
  geom_point(aes(shape = excludes_0), color = "black", size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.147, 0.147),
             linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = c(-0.33, 0.33),
             linetype = "dotted", color = "grey60") +
  scale_fill_manual(
    values = c(negligible = "grey80", small = "#fdcc8a",
               medium = "#fc8d59", large = "#d7301f"),
    drop = FALSE
  ) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "CI excludes 0", `FALSE` = "CI crosses 0"),
                     name = NULL) +
  labs(
    title = "Cliff's delta: Mesenchymal vs Non-Mesenchymal (ALT-high only)",
    subtitle = sprintf("Stratified to ALT-high subset (n=%d mes vs %d non-mes). Bars: 95%% bootstrap CI.",
                       sum(df_alt_high$cancer_group == grp_mes),
                       sum(df_alt_high$cancer_group == grp_non)),
    x = "Cliff's delta", y = NULL, fill = "Effect size"
  ) +
  theme_bw()

plt_effect_ranking
# ---- 3. Distribution plots for top-k features ---------------

top_features <- group_diff_numeric %>%
  filter(cliff_mag %in% c("small", "medium", "large")) %>%
  pull(feature)

if (length(top_features) == 0) {
  message("No features with non-negligible effect size. Showing top 6 by |delta|.")
  top_features <- group_diff_numeric %>%
    slice_max(abs(cliff_delta), n = 6) %>%
    pull(feature)
}

dist_plot <- function(v) {
  df_group %>%
    ggplot(aes(x = cancer_group, y = .data[[v]], fill = cancer_group)) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = cancer_group),
                width = 0.08, alpha = 0.5, size = 0.8) +
    scale_fill_manual(values = grp_colors) +
    scale_color_manual(values = grp_colors_dark) +
    labs(x = NULL, y = v, title = v) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.text.x = element_text(size = 8))
}

dist_plots <- map(top_features, dist_plot)
wrap_plots(dist_plots, ncol = 3) +
  plot_annotation(
    title    = "Top features by effect size — distribution by cancer group",
    subtitle = sprintf("%s in red. Note jitter density: %d vs %d.",
                       grp_mes,
                       sum(df_group$cancer_group == grp_mes),
                       sum(df_group$cancer_group == grp_non))
  )


# ---- 4. ALT-status composition by group ---------------------

alt_table <- df_group %>%
  count(cancer_group, alt_status) %>%
  group_by(cancer_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

plt_alt <- alt_table %>%
  ggplot(aes(x = cancer_group, y = prop, fill = alt_status)) +
  geom_col(position = "stack", alpha = 0.85) +
  geom_text(aes(label = paste0(n, "\n(", scales::percent(prop, 1), ")")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("ALT-low" = "grey70", "ALT-high" = "firebrick")) +
  labs(title = "ALT status composition by cancer group",
       x = NULL, y = "Proportion", fill = NULL) +
  theme_bw()

plt_alt

# Stratified bootstrap CI on the difference in P(ALT-high)
set.seed(21)
alt_high_mes <- as.numeric(df_group$alt_status[df_group$cancer_group == grp_mes] == "ALT-high")
alt_high_non <- as.numeric(df_group$alt_status[df_group$cancer_group == grp_non] == "ALT-high")

boot_diff <- replicate(2000, {
  s_mes <- sample(alt_high_mes, length(alt_high_mes), replace = TRUE)
  s_non <- sample(alt_high_non, length(alt_high_non), replace = TRUE)
  mean(s_mes) - mean(s_non)
})
boot_diff <- boot_diff[is.finite(boot_diff)]

p_mes_obs <- mean(alt_high_mes)
p_non_obs <- mean(alt_high_non)

cat("\nALT-high prevalence by cancer group\n")
cat(sprintf("  %s: %.1f%% (%d/%d)\n", grp_mes,
            100 * p_mes_obs, sum(alt_high_mes), length(alt_high_mes)))
cat(sprintf("  %s: %.1f%% (%d/%d)\n", grp_non,
            100 * p_non_obs, sum(alt_high_non), length(alt_high_non)))
cat(sprintf("\nDifference (mes - non):\n"))
cat(sprintf("  observed: %+.3f\n", p_mes_obs - p_non_obs))
cat(sprintf("  95%% stratified bootstrap CI: [%+.3f, %+.3f]\n",
            quantile(boot_diff, 0.025), quantile(boot_diff, 0.975)))


# ---- 5. Response variable by group --------------------------

plt_response <- df_group %>%
  ggplot(aes(x = cancer_group, y = tf_primary_rate, fill = cancer_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.9) +
  scale_fill_manual(values = grp_colors) +
  labs(title = "tf_primary_rate (log1p) by cancer group",
       subtitle = "Is the prediction target itself shifted?",
       x = NULL, y = "tf_primary_rate (log1p)") +
  theme_bw() +
  theme(legend.position = "none")

plt_response_alt <- df_group %>%
  ggplot(aes(x = cancer_group, y = tf_primary_rate, fill = cancer_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  facet_wrap(~ alt_status) +
  scale_fill_manual(values = grp_colors) +
  labs(title = "tf_primary_rate by group, stratified by ALT status",
       subtitle = "Check whether group effect survives ALT stratification",
       x = NULL, y = "tf_primary_rate (log1p)") +
  theme_bw() +
  theme(legend.position = "none")

plt_response / plt_response_alt


plt_response <- df_group %>%
  ggplot(aes(x = cancer_group, y = tf_blood_rate, fill = cancer_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.9) +
  scale_fill_manual(values = grp_colors) +
  labs(title = "tf_blood_rate (log1p) by cancer group",
       subtitle = "Is the prediction target itself shifted?",
       x = NULL, y = "tf_blood_rate (log1p)") +
  theme_bw() +
  theme(legend.position = "none")

plt_response_alt <- df_group %>%
  ggplot(aes(x = cancer_group, y = tf_blood_rate, fill = cancer_group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  facet_wrap(~ alt_status) +
  scale_fill_manual(values = grp_colors) +
  labs(title = "tf_blood_rate by group, stratified by ALT status",
       subtitle = "Check whether group effect survives ALT stratification",
       x = NULL, y = "tf_blood_rate (log1p)") +
  theme_bw() +
  theme(legend.position = "none")

plt_response / plt_response_alt


# ---- 6. Bivariate view -------------------------------------

top2 <- group_diff_plot$feature[order(-abs(group_diff_plot$cliff_delta))][1:2]
top2 <- as.character(top2)

plt_bivariate <- df_group %>%
  ggplot(aes(x = .data[[top2[1]]], y = .data[[top2[2]]],
             color = cancer_group, size = cancer_group,
             alpha = cancer_group)) +
  geom_point() +
  scale_color_manual(values = grp_colors) +
  scale_size_manual(values = grp_sizes) +
  scale_alpha_manual(values = grp_alphas) +
  labs(title = paste0("Bivariate distribution: ", top2[1], " vs ", top2[2]),
       subtitle = "Do mesenchymals cluster in a specific region of the feature space?",
       color = NULL, size = NULL, alpha = NULL) +
  theme_bw()

plt_bivariate


# ---- 7. Correlation structure differences -------------------

cor_mes <- df_group %>%
  filter(cancer_group == grp_mes) %>%
  select(all_of(numeric_vars)) %>%
  cor(use = "pairwise.complete.obs")

cor_non <- df_group %>%
  filter(cancer_group == grp_non) %>%
  select(all_of(numeric_vars)) %>%
  cor(use = "pairwise.complete.obs")

cor_diff <- cor_mes - cor_non

cor_diff_long <- cor_diff %>%
  as.data.frame() %>%
  rownames_to_column("feat1") %>%
  pivot_longer(-feat1, names_to = "feat2", values_to = "diff") %>%
  filter(feat1 < feat2) %>%
  arrange(desc(abs(diff)))

cat("\nTop pairwise correlation differences (mes - non):\n")
print(head(cor_diff_long, 15))

plt_cor_diff <- cor_diff %>%
  as.data.frame() %>%
  rownames_to_column("feat1") %>%
  pivot_longer(-feat1, names_to = "feat2", values_to = "diff") %>%
  ggplot(aes(x = feat1, y = feat2, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation structure difference",
       subtitle = sprintf("cor(%s) - cor(%s). Red = stronger in mesenchymal.",
                          grp_mes, grp_non),
       x = NULL, y = NULL, fill = "Δ cor") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))

plt_cor_diff