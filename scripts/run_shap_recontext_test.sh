#!/bin/bash
#PBS -o /ngs/iflores/andrea/logs/shap_recon_${PBS_JOBID}.log
#PBS -e /ngs/iflores/andrea/logs/shap_recon_${PBS_JOBID}.error

# Recontextualization explainer, TEST on 5 rows only. Logit-space output.
# Goal: measure GPU wall-time and read the diagnostics to confirm explained class.
# Submit with GPU:  qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_recontext_test.sh

source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd /ngs/iflores/andrea
mkdir -p logs output_recontext_test

python SHAP_tabpfn_local.py \
  --explainer  recontextualization \
  --input-dir  /ngs/iflores/andrea/input \
  --output-dir /ngs/iflores/andrea/output_recontext_test \
  --budget     128 \
  --n-explain  5 \
  --n-jobs     16 \
  --class-index 1 \
  --seed       42
