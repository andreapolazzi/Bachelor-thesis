#!/bin/bash
#
# PRODUCTION SHAP run on the standalone 'carmela' box (RTX 3080 GPU, /ngs mounted).
# Computes BOTH explainers for ALL test rows:
#   STEP 1: imputation          -> output_imputation_full/  (probability space, self-verified)
#   STEP 2: recontextualization -> output_recontext_full/   (TabPFN-native logit-like space)
# Both are verified-correct for class 1 (ALT-high); see diag_recontext.py.
#
# Every line of stdout+stderr is tee'd to a timestamped log under logs/, and each
# step prints its own exit code, so any failure is fully traceable.
#
# Run directly (no scheduler). Survive SSH disconnects with nohup:
#     cd /ngs/iflores/andrea
#     nohup bash run_shap_full.sh > /tmp/shap_full.out 2>&1 &
#     tail -f /tmp/shap_full.out

set -u

BASE=/ngs/iflores/andrea
mkdir -p "$BASE/logs" "$BASE/output_imputation_full" "$BASE/output_recontext_full"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_full_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP full (both explainers)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

source /ngs/software/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate /ngs/software/conda/envs/tabpfn2 \
  || micromamba activate /ngs/software/conda/envs/tabpfn2 \
  || { echo "ERROR: could not activate tabpfn2 env"; exit 1; }
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

run_step "imputation-full" \
  --explainer   imputation \
  --input-dir   "$BASE/input" \
  --output-dir  "$BASE/output_imputation_full" \
  --budget      256 \
  --n-jobs      4 \
  --class-index 1 \
  --seed        42
IMP_RC=$?

run_step "recontext-full" \
  --explainer   recontextualization \
  --input-dir   "$BASE/input" \
  --output-dir  "$BASE/output_recontext_full" \
  --budget      256 \
  --n-jobs      4 \
  --class-index 1 \
  --seed        42
RECON_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY"
echo "#   imputation-full exit code : $IMP_RC   (0 = OK)"
echo "#   recontext-full  exit code : $RECON_RC (0 = OK)"
echo "#   finished $(date)"
echo "#   full log: $LOG"
echo "######################################################################"
