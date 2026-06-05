#!/bin/bash
#
# PRODUCTION OOF SHAP run on the chiron HPC GPU queue (16 cores).
# Both explainers, 5-fold out-of-fold, all samples, budget 256.
#   imputation     -> output_imputation_full2/   (probability space)
#   recontextualiz -> output_recontext_full2/    (logit-like space)
#
# Logging is internal (tee to a timestamped file) -- do NOT rely on #PBS -o with
# ${PBS_JOBID} (it does not expand in directive lines and the log is lost).
#
# Submit (GPU, 16 cores):
#   qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_full2.sh

set -u

BASE=/ngs/iflores/andrea
mkdir -p "$BASE/logs" "$BASE/output_imputation_full2" "$BASE/output_recontext_full2"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_full2_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP OOF full2 (both explainers)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# env activation MUST precede python (else ModuleNotFoundError: tabpfn)
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

echo "host        : $(hostname)"
echo "which python: $(which python)"
nvidia-smi || echo "(no nvidia-smi / no GPU visible)"
echo "-----------------------------------"

run_step () {
  local name="$1"; shift
  echo
  echo "=================================================================="
  echo "START $name  |  $(date)"
  echo "cmd: python -u SHAP_tabpfn_local.py $*"
  echo "=================================================================="
  python -u SHAP_tabpfn_local.py "$@"
  local rc=$?
  echo "------------------------------------------------------------------"
  echo "END   $name  |  exit code = $rc  |  $(date)"
  echo "=================================================================="
  return $rc
}

run_step "imputation-oof" \
  --explainer   imputation \
  --input-dir   "$BASE/shap_input_cv" \
  --output-dir  "$BASE/output_imputation_full2" \
  --folds       "$BASE/shap_input_cv/folds.csv" \
  --budget      256 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
IMP_RC=$?

run_step "recontext-oof" \
  --explainer   recontextualization \
  --input-dir   "$BASE/shap_input_cv" \
  --output-dir  "$BASE/output_recontext_full2" \
  --folds       "$BASE/shap_input_cv/folds.csv" \
  --budget      256 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
RECON_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY"
echo "#   imputation-oof exit code : $IMP_RC   (0 = OK)"
echo "#   recontext-oof  exit code : $RECON_RC (0 = OK)"
echo "#   finished $(date)  |  log: $LOG"
echo "######################################################################"
