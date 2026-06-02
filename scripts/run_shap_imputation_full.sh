#!/bin/bash
#PBS -o /ngs/iflores/andrea/logs/shap_imp_${PBS_JOBID}.log
#PBS -e /ngs/iflores/andrea/logs/shap_imp_${PBS_JOBID}.error

# Imputation explainer, FULL test set. Probability-space, self-verified output.
# Submit with GPU:  qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_imputation_full.sh

source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd /ngs/iflores/andrea
mkdir -p logs output_imputation_full

python SHAP_tabpfn_local.py \
  --explainer  imputation \
  --input-dir  /ngs/iflores/andrea/input \
  --output-dir /ngs/iflores/andrea/output_imputation_full \
  --budget     256 \
  --n-jobs     16 \
  --class-index 1 \
  --seed       42
