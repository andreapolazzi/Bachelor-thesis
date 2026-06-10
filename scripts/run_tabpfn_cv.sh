#!/bin/bash
#
# TabPFN cross-validation (classification + both TF-rate regressions) on the chiron
# HPC GPU queue. Local `tabpfn` (no cloud API), shared 5x5 stratified folds, out-of-fold
# predictions + per-fold metrics — so the TabPFN rows are directly comparable to the
# in-R models (logistic reg / RF / LASSO for classification; LM / LASSO / RF for
# regression) scored on the SAME folds in notebooks 20 and 21.
#
# Three jobs:
#   classification  shap_input_cmp_cv -> tabpfn_cv_oof.csv      (+ _summary.csv)
#   tTF regression  cmp_reg_ttf_cv    -> tabpfn_reg_ttf_oof.csv (+ _summary.csv)
#   bTF regression  cmp_reg_btf_cv    -> tabpfn_reg_btf_oof.csv (+ _summary.csv)
#
# Each input dir holds X_full.csv, y_full.csv, folds_5x5.csv, produced locally by
#   Rscript scripts/export_cv_inputs.R          # (or by rendering notebooks 20 & 21)
# and synced into $BASE alongside tabpfn_cv_oof.py. After the run, sync the *_oof.csv
# back to the project outputs/ and re-render notebooks 20 & 21 for the final tables.
#
# Submit (GPU, 4 cores — n=16 OOM'd the RTX 3080 in the SHAP runs):
#   qsubmit.pl -g 1 -n 4 -s /ngs/iflores/andrea/run_tabpfn_cv.sh
# or on a no-scheduler box (carmela):
#   cd "$BASE" && nohup bash run_tabpfn_cv.sh >/dev/null 2>&1 &

set -u

BASE="${TABPFN_BASE:-/ngs/iflores/andrea}"
mkdir -p "$BASE/logs"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/tabpfn_cv_${TS}.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# TabPFN CV (classification + tTF + bTF)  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# env activation MUST precede python (else ModuleNotFoundError: tabpfn)
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

# reduce GPU memory fragmentation (same as the SHAP re-runs)
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "host        : $(hostname)"
echo "which python: $(which python)"
echo "PYTORCH_CUDA_ALLOC_CONF: $PYTORCH_CUDA_ALLOC_CONF"
nvidia-smi || echo "(no nvidia-smi / no GPU visible)"
echo "-----------------------------------"

# pre-flight: the export dirs (X_full / y_full / folds_5x5) must be present
missing=0
for d in shap_input_cmp_cv cmp_reg_ttf_cv cmp_reg_btf_cv; do
  for f in X_full.csv y_full.csv folds_5x5.csv; do
    [ -f "$BASE/$d/$f" ] || { echo "MISSING: $BASE/$d/$f"; missing=1; }
  done
done
if [ "$missing" -ne 0 ]; then
  echo
  echo "Input dirs incomplete. Generate them locally and sync into \$BASE:"
  echo "  Rscript scripts/export_cv_inputs.R"
  echo "  rsync -a outputs/shap_input_cmp_cv outputs/cmp_reg_ttf_cv outputs/cmp_reg_btf_cv  $BASE/"
  exit 2
fi

run_step () {
  local name="$1"; shift
  echo
  echo "=================================================================="
  echo "START $name  |  $(date)"
  echo "cmd: python -u tabpfn_cv_oof.py $*"
  echo "=================================================================="
  python -u tabpfn_cv_oof.py "$@"
  local rc=$?
  echo "------------------------------------------------------------------"
  echo "END   $name  |  exit code = $rc  |  $(date)"
  echo "=================================================================="
  return $rc
}

run_step "classification" \
  --task        classification \
  --input-dir   "$BASE/shap_input_cmp_cv" \
  --output      "$BASE/tabpfn_cv_oof.csv" \
  --n-jobs      4 \
  --seed        21
CLF_RC=$?

run_step "ttf-regression" \
  --task        regression \
  --input-dir   "$BASE/cmp_reg_ttf_cv" \
  --output      "$BASE/tabpfn_reg_ttf_oof.csv" \
  --n-jobs      4 \
  --seed        21
TTF_RC=$?

run_step "btf-regression" \
  --task        regression \
  --input-dir   "$BASE/cmp_reg_btf_cv" \
  --output      "$BASE/tabpfn_reg_btf_oof.csv" \
  --n-jobs      4 \
  --seed        21
BTF_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY  (0 = OK)"
echo "#   classification  exit code : $CLF_RC"
echo "#   ttf-regression  exit code : $TTF_RC"
echo "#   btf-regression  exit code : $BTF_RC"
echo "#"
echo "#   outputs: tabpfn_cv_oof.csv  tabpfn_reg_ttf_oof.csv  tabpfn_reg_btf_oof.csv"
echo "#            (+ matching *_summary.csv with mean/SD of every metric)"
echo "#   next: sync the *_oof.csv back to the project outputs/ and re-render"
echo "#         notebooks 20 (classification) and 21 (regression)."
echo "#   finished $(date)  |  log: $LOG"
echo "######################################################################"
