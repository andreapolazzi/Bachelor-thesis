# Graph Report - .  (2026-06-11)

## Corpus Check
- 54 files · ~57,183 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 132 nodes · 188 edges · 18 communities (11 shown, 7 thin omitted)
- Extraction: 90% EXTRACTED · 10% INFERRED · 0% AMBIGUOUS · INFERRED: 18 edges (avg confidence: 0.85)
- Token cost: 265,118 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_OOF SHAP CV Harness|OOF SHAP CV Harness]]
- [[_COMMUNITY_TabPFN Pipeline Scripts|TabPFN Pipeline Scripts]]
- [[_COMMUNITY_Local SHAP Explainer|Local SHAP Explainer]]
- [[_COMMUNITY_SHAP Importance Bar Figure|SHAP Importance Bar Figure]]
- [[_COMMUNITY_ALT Status & Confounding|ALT Status & Confounding]]
- [[_COMMUNITY_TF-Rate Regression|TF-Rate Regression]]
- [[_COMMUNITY_SHAP Beeswarm Figure|SHAP Beeswarm Figure]]
- [[_COMMUNITY_TERT & Feature Selection|TERT & Feature Selection]]
- [[_COMMUNITY_SHAP Dependence Figure|SHAP Dependence Figure]]
- [[_COMMUNITY_Regression SHAP Runner|Regression SHAP Runner]]
- [[_COMMUNITY_TabPFN CV Runner|TabPFN CV Runner]]
- [[_COMMUNITY_SHAP Full Runner|SHAP Full Runner]]
- [[_COMMUNITY_SHAP Full Runner 2|SHAP Full Runner 2]]
- [[_COMMUNITY_SHAP Imputation Runner|SHAP Imputation Runner]]
- [[_COMMUNITY_SHAP Test Runner|SHAP Test Runner]]
- [[_COMMUNITY_SHAP Both Runner|SHAP Both Runner]]
- [[_COMMUNITY_SHAP Imputation Full Runner|SHAP Imputation Full Runner]]
- [[_COMMUNITY_SHAP Recontext Runner|SHAP Recontext Runner]]

## God Nodes (most connected - your core abstractions)
1. `10_TabPFN (TabPFN classification + regression OOF SHAP)` - 13 edges
2. `TabPFN aggregate SHAP feature importance bar chart` - 11 edges
3. `PCAWG_primary.xlsx (primary tumor dataset)` - 10 edges
4. `TabPFN SHAP beeswarm summary plot` - 9 edges
5. `main()` - 8 edges
6. `run_oof()` - 7 edges
7. `13_1_tTF_regression_+blood (tTF regression incl. blood TF)` - 7 edges
8. `15_1_bTF_regression_+tumor (bTF regression incl. tumor TF)` - 6 edges
9. `20_classification_comparison (fair ALT classifier comparison)` - 6 edges
10. `21_regression_comparison (fair TF-rate regressor comparison)` - 6 edges

## Surprising Connections (you probably didn't know these)
- `10_TabPFN (TabPFN classification + regression OOF SHAP)` --references--> `carmela GPU box (local TabPFN/SHAP, no scheduler)`  [EXTRACTED]
  analysis/primary/10_TabPFN.qmd → docs/superpowers/specs/2026-06-08-regression-oof-shap-design.md
- `10_TabPFN (TabPFN classification + regression OOF SHAP)` --implements--> `Out-of-fold SHAP (each sample explained by model not trained on it)`  [EXTRACTED]
  analysis/primary/10_TabPFN.qmd → docs/superpowers/specs/2026-06-05-blood-inclusive-oof-shap-design.md
- `main()` --calls--> `Path`  [INFERRED]
  scripts/SHAP_tabpfn_local.py → scripts/test_oof_smoke.py
- `main()` --calls--> `Path`  [INFERRED]
  scripts/run_tabpfn_shap.py → scripts/test_oof_smoke.py
- `mesenchymal_groups (Mesenchymal vs Non-Mesenchymal comparison)` --semantically_similar_to--> `6_TVRs_comparison (TVR singleton distance primary vs metastatic)`  [INFERRED] [semantically similar]
  analysis/primary/mesenchymal_groups.qmd → analysis/metastatic/6_TVRs_comparison.qmd

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Shared 5x5 stratified CV harness across classification, regression, and SHAP** — concept_make_folds, primary_20_classification_comparison, primary_21_regression_comparison, primary_10_tabpfn [EXTRACTED 0.95]
- **OOF SHAP pipeline (spec to plan to notebook to local script)** — specs_2026_06_05_blood_inclusive_oof_shap_design, plans_2026_06_05_blood_inclusive_oof_shap, primary_10_tabpfn, concept_shap_tabpfn_local_py, concept_reconstruct_shapviz [EXTRACTED 0.90]
- **Two-part TF-rate regression notebooks (tTF/bTF, with/without TERT and blood/tumor)** — primary_13_ttf_regression_primary, primary_13_1_ttf_regression_blood, primary_15_btf_regression_primary, primary_15_1_btf_regression_tumor, primary_16_ttf_regression_notert, primary_17_btf_regression_notert [INFERRED 0.85]

## Communities (18 total, 7 thin omitted)

### Community 0 - "OOF SHAP CV Harness"
Cohesion: 0.23
Nodes (16): carmela GPU box (local TabPFN/SHAP, no scheduler), scripts/clf_metrics.R (fold_metrics / avg_precision), features14 (14-feature telomere set incl. tf_blood_rate), Imputation + recontextualization SHAP explainers, scripts/make_folds.R (shared 5x5 stratified CV folds, seed 21), Out-of-fold SHAP (each sample explained by model not trained on it), scripts/reconstruct_shapviz.R (load_shap / load_meta), scripts/reg_metrics.R (reg_fold_metrics: RMSE/MAE/Spearman) (+8 more)

### Community 1 - "TabPFN Pipeline Scripts"
Cohesion: 0.19
Nodes (12): Path, main(), sig(), main(), main(), make_model(), Per-fold + summary metrics, matching clf_metrics.R / reg_metrics.R so the     Ta, score() (+4 more)

### Community 2 - "Local SHAP Explainer"
Cohesion: 0.25
Nodes (13): build_explainer(), extract_shap(), fit_model(), load_data(), main(), make_predict_fn(), parse_args(), Fit and return a TabPFN model for the requested task. (+5 more)

### Community 3 - "SHAP Importance Bar Figure"
Cohesion: 0.18
Nodes (13): CTAGGG_singleton_dist, TabPFN aggregate SHAP feature importance bar chart, GTAGGG_singleton_dist, mean absolute SHAP value, telomeric variant repeat singleton_dist motifs, TAAGGG_singleton_dist, TabPFN classification model, telomere_content_log2 (+5 more)

### Community 4 - "ALT Status & Confounding"
Cohesion: 0.23
Nodes (12): alt_status (ALT-low / ALT-high label), cancer_group (Mesenchymal vs Non-Mesenchymal origin), Cancer-type confounding of ALT/feature signal, metastatic_red_edit_mesenchymal.xlsx, metastatic_red_edit_singleton_dist.xlsx (long-format metastatic), PCAWG_primary+genes.xlsx (primary dataset with gene columns), TERT_FPKM spurious positive predictor (shared variance + missingness artifact), TVR singleton_dist motifs (9 telomere variant repeats) (+4 more)

### Community 5 - "TF-Rate Regression"
Cohesion: 0.35
Nodes (12): log1p transform of skewed zero-inflated rates, PCAWG_primary.xlsx (primary tumor dataset), tf_blood_rate / bTF (blood telomere fusion rate), tf_primary_rate / tTF (tumor telomere fusion rate), Two-part model (presence logistic + intensity Gamma GLM, Duan smearing), 13_1_tTF_regression_+blood (tTF regression incl. blood TF), 13_tTF_regression_primary (tTF regression), 15_1_bTF_regression_+tumor (bTF regression incl. tumor TF) (+4 more)

### Community 6 - "SHAP Beeswarm Figure"
Cohesion: 0.22
Nodes (11): TabPFN SHAP beeswarm summary plot, SHAP value (impact on model output), TVR singleton distribution feature family, TabPFN model, telomere_content_log2, telomere_insertion_rate, TERT_FPKM, tf_rate (+3 more)

### Community 7 - "TERT & Feature Selection"
Cohesion: 0.43
Nodes (7): scripts/clean_columns.R (clean_pcawg_data), PCAWG_variables_PCA_woK.xlsx (raw PCA variables), TERT_FPKM (TERT expression), 11_Boruta_feat_imp (Boruta feature importance), 12_dataset_balancing (ROSE/SMOTE/class-weight comparison), 14_TERT_missingness (TERT_FPKM missingness analysis), 16_tTF_regression_noTERT (tTF regression without TERT)

### Community 8 - "SHAP Dependence Figure"
Cohesion: 0.50
Nodes (5): TabPFN SHAP dependence plot for telomere_insertion_rate, SHAP value, TabPFN model, telomere_insertion_rate, tf_rate

### Community 9 - "Regression SHAP Runner"
Cohesion: 0.67
Nodes (3): PYTORCH_CUDA_ALLOC_CONF, run_step(), run_shap_reg.sh script

### Community 10 - "TabPFN CV Runner"
Cohesion: 0.67
Nodes (3): PYTORCH_CUDA_ALLOC_CONF, run_step(), run_tabpfn_cv.sh script

## Knowledge Gaps
- **23 isolated node(s):** `run_shap_both.sh script`, `run_shap_imp_full2.sh script`, `PYTORCH_CUDA_ALLOC_CONF`, `run_shap_imputation_full.sh script`, `run_shap_recontext_test.sh script` (+18 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **7 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `PCAWG_primary.xlsx (primary tumor dataset)` connect `TF-Rate Regression` to `OOF SHAP CV Harness`, `ALT Status & Confounding`?**
  _High betweenness centrality (0.055) - this node is a cross-community bridge._
- **Why does `10_TabPFN (TabPFN classification + regression OOF SHAP)` connect `OOF SHAP CV Harness` to `TF-Rate Regression`, `TERT & Feature Selection`?**
  _High betweenness centrality (0.048) - this node is a cross-community bridge._
- **Why does `Path` connect `TabPFN Pipeline Scripts` to `Local SHAP Explainer`?**
  _High betweenness centrality (0.034) - this node is a cross-community bridge._
- **What connects `Fit and return a TabPFN model for the requested task.`, `Return (predict_fn, proba_col).      For classification, predict_fn returns pred`, `Return (explainer, space) for the chosen paradigm, with X_bg as background.` to the rest of the system?**
  _33 weakly-connected nodes found - possible documentation gaps or missing edges._