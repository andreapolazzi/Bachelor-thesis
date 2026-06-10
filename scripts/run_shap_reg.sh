#!/bin/bash
#
# OOF SHAP for the TabPFN REGRESSION tasks (tTF = tf_primary_rate, bTF = tf_blood_rate).
# Both explainers, single stratified 5-fold out-of-fold, all samples, budget 256.
#   imputation -> *_imp/   (target log1p space, additive & self-verified)
#   recontext  -> *_rec/   (recontextualized target space, author-recommended cross-check)
#
# TabPFN is GPU-bound: CPU parallelism gives no speedup and just multiplies GPU memory,
# so we use --n-jobs 4 with expandable CUDA segments (imputation OOM'd at 16 in full2).
#
# Designed for the carmela standalone GPU box (no scheduler), launched with nohup:
#   cd "$SHAP_BASE"
#   nohup bash run_shap_reg.sh >/dev/null 2>&1 &
# It also works under PBS on chiron:
#   qsubmit.pl -g 1 -n 4 -s /ngs/iflores/andrea/run_shap_reg.sh
#
# Set SHAP_BASE to the dir that holds SHAP_tabpfn_local.py + the shap_input_*_cv/ folders.

set -u

BASE="${SHAP_BASE:-/ngs/iflores/andrea}"
mkdir -p "$BASE/logs" \
  "$BASE/shap_output_ttf_imp" "$BASE/shap_output_ttf_rec" \
  "$BASE/shap_output_btf_imp" "$BASE/shap_output_btf_rec"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_reg_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP OOF regression (tTF & bTF, both explainers)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# env activation MUST precede python (else ModuleNotFoundError: tabpfn).
# Adjust this block for carmela if its env path differs.
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

# reduce GPU memory fragmentation (avoids the OOM seen in full2)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "host        : $(hostname)"
echo "which python: $(which python)"
echo "PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"
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

# --- tTF (tf_primary_rate) -------------------------------------------------
run_step "ttf-imputation-oof" \
  --task        regression \
  --explainer   imputation \
  --input-dir   "$BASE/shap_input_ttf_cv" \
  --output-dir  "$BASE/shap_output_ttf_imp" \
  --folds       "$BASE/shap_input_ttf_cv/folds.csv" \
  --budget      256 \
  --n-jobs      4 \
  --seed        42
TTF_IMP_RC=$?

run_step "ttf-recontext-oof" \
  --task        regression \
  --explainer   recontextualization \
  --input-dir   "$BASE/shap_input_ttf_cv" \
  --output-dir  "$BASE/shap_output_ttf_rec" \
  --folds       "$BASE/shap_input_ttf_cv/folds.csv" \
  --budget      256 \
  --n-jobs      4 \
  --seed        42
TTF_REC_RC=$?

# --- bTF (tf_blood_rate) ---------------------------------------------------
run_step "btf-imputation-oof" \
  --task        regression \
  --explainer   imputation \
  --input-dir   "$BASE/shap_input_btf_cv" \
  --output-dir  "$BASE/shap_output_btf_imp" \
  --folds       "$BASE/shap_input_btf_cv/folds.csv" \
  --budget      256 \
  --n-jobs      4 \
  --seed        42
BTF_IMP_RC=$?

run_step "btf-recontext-oof" \
  --task        regression \
  --explainer   recontextualization \
  --input-dir   "$BASE/shap_input_btf_cv" \
  --output-dir  "$BASE/shap_output_btf_rec" \
  --folds       "$BASE/shap_input_btf_cv/folds.csv" \
  --budget      256 \
  --n-jobs      4 \
  --seed        42
BTF_REC_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY  (0 = OK)"
echo "#   ttf-imputation-oof exit code : $TTF_IMP_RC"
echo "#   ttf-recontext-oof  exit code : $TTF_REC_RC"
echo "#   btf-imputation-oof exit code : $BTF_IMP_RC"
echo "#   btf-recontext-oof  exit code : $BTF_REC_RC"
echo "#   finished $(date)  |  log: $LOG"
echo "######################################################################"
