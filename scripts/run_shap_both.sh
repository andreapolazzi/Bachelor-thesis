#!/bin/bash
#PBS -o /ngs/iflores/andrea/logs/shap_both_${PBS_JOBID}.log
#PBS -e /ngs/iflores/andrea/logs/shap_both_${PBS_JOBID}.error

# Runs BOTH explainers in a single job (queue once, hold the GPU for both):
#   1) recontextualization TEST on 5 rows   -> output_recontext_test/   (logit space)
#   2) imputation FULL on all rows          -> output_imputation_full/  (probability space)
# The recontext test runs first because it is short; if it crashes, the full run
# still proceeds (we do not abort the job on its failure).
# Submit with GPU:  qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_both.sh

source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd /ngs/iflores/andrea
mkdir -p logs output_recontext_test output_imputation_full

echo "=================================================================="
echo "STEP 1/2: recontextualization test (5 rows)"
echo "=================================================================="
python SHAP_tabpfn_local.py \
  --explainer  recontextualization \
  --input-dir  /ngs/iflores/andrea/input \
  --output-dir /ngs/iflores/andrea/output_recontext_test \
  --budget     128 \
  --n-explain  5 \
  --n-jobs     16 \
  --class-index 1 \
  --seed       42 \
  || echo "WARNING: recontextualization step failed; continuing to full imputation run."

echo
echo "=================================================================="
echo "STEP 2/2: imputation full run (all rows)"
echo "=================================================================="
python SHAP_tabpfn_local.py \
  --explainer  imputation \
  --input-dir  /ngs/iflores/andrea/input \
  --output-dir /ngs/iflores/andrea/output_imputation_full \
  --budget     256 \
  --n-jobs     16 \
  --class-index 1 \
  --seed       42

echo
echo "ALL DONE."
