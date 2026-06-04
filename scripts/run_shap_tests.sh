#!/bin/bash
#PBS -N shap_tests
#
# Two SMALL, fully-traceable SHAP test runs in one job.
#   STEP 1: imputation (probability space, self-verified)  -> output_imp_test/
#   STEP 2: recontextualization (logit space, diagnostic)  -> output_recon_test/
#
# Logging is done INSIDE the script to a timestamped file, NOT via #PBS -o/-e.
# Reason: ${PBS_JOBID} does not expand in #PBS directive lines, so those logs
# get lost. Here every line of stdout+stderr is tee'd to:
#     /ngs/iflores/andrea/logs/shap_tests_<timestamp>.log
# so a crash is always traceable.
#
# Submit with GPU:  qsubmit.pl -g 1 -n 16 -s /ngs/iflores/andrea/run_shap_tests.sh

set -u

BASE=/ngs/iflores/andrea
mkdir -p "$BASE/logs" "$BASE/output_imp_test" "$BASE/output_recon_test"

TS=$(date +%Y%m%d_%H%M%S)
LOG="$BASE/logs/shap_tests_${TS}.log"

# Send EVERYTHING (stdout + stderr) to the log file AND to the console.
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# SHAP tests  |  started $(date)"
echo "# log file: $LOG"
echo "######################################################################"

# --- environment ---
source /ngs/software/micromamba/env.sh
micromamba activate /ngs/software/conda/envs/tabpfn2
cd "$BASE"

echo "host        : $(hostname)"
echo "which python: $(which python)"
echo "PBS_JOBID   : ${PBS_JOBID:-<not set>}"
echo "--- GPU visibility (nvidia-smi) ---"
nvidia-smi || echo "(no nvidia-smi / no GPU visible)"
echo "-----------------------------------"

# Helper: run one step, always print its exit code, never abort the whole job.
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

# STEP 1: imputation, small (10 rows). This is the trustworthy, self-verified one.
run_step "imputation-test" \
  --explainer   imputation \
  --input-dir   "$BASE/input" \
  --output-dir  "$BASE/output_imp_test" \
  --budget      128 \
  --n-explain   10 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
IMP_RC=$?

# STEP 2: recontextualization, small (10 rows). Diagnostic only (logit space).
run_step "recontext-test" \
  --explainer   recontextualization \
  --input-dir   "$BASE/input" \
  --output-dir  "$BASE/output_recon_test" \
  --budget      128 \
  --n-explain   10 \
  --n-jobs      16 \
  --class-index 1 \
  --seed        42
RECON_RC=$?

echo
echo "######################################################################"
echo "# SUMMARY"
echo "#   imputation-test exit code : $IMP_RC   (0 = OK)"
echo "#   recontext-test  exit code : $RECON_RC (0 = OK)"
echo "#   finished $(date)"
echo "#   full log: $LOG"
echo "######################################################################"
