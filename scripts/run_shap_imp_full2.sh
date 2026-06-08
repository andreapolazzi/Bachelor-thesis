#!/bin/bash
#
# RE-RUN of the imputation OOF SHAP only (recontext already succeeded).
# The full2 run OOM'd on the RTX 3080 because --n-jobs 16 spawned 16 TabPFN
# GPU model copies that saturated the 10 GB card. For GPU-bound TabPFN, CPU
# parallelism gives no speedup and just multiplies GPU memory, so we drop to
# --n-jobs 4 and enable expandable CUDA segments to reduce fragmentation.
#
#   imputation -> output_imputation_full2/   (probability space)
#
# Submit (GPU, fewer cores fine since GPU-bound):
#   qsubmit.pl -g 1 -n 4 -s /ngs/iflores/andrea/run_shap_imp_full2.sh

set -u

BASE=/ngs/iflores/andrea
mkdir -p "$BASE/logs" "$BASE/output_imputation_full2"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_imp_full2_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP OOF imputation re-run (n-jobs 4)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# env activation MUST precede python (else ModuleNotFoundError: tabpfn)
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

# reduce GPU memory fragmentation (suggested by the OOM error)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "host        : $(hostname)"
echo "which python: $(which python)"
echo "PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"
nvidia-smi || echo "(no nvidia-smi / no GPU visible)"
echo "-----------------------------------"

echo
echo "=================================================================="
echo "START imputation-oof (n-jobs 4)  |  $(date)"
echo "=================================================================="
python -u SHAP_tabpfn_local.py \
  --explainer   imputation \
  --input-dir   "$BASE/shap_input_cv" \
  --output-dir  "$BASE/output_imputation_full2" \
  --folds       "$BASE/shap_input_cv/folds.csv" \
  --budget      256 \
  --n-jobs      4 \
  --class-index 1 \
  --seed        42
RC=$?
echo "------------------------------------------------------------------"
echo "END   imputation-oof  |  exit code = $RC  |  $(date)"
echo "######################################################################"
echo "#   imputation-oof exit code : $RC   (0 = OK)  |  log: $LOG"
echo "######################################################################"
