#!/bin/bash
#
# Two SMALL, fully-traceable SHAP test runs, for the STANDALONE 'carmela' box
# (no PBS/Slurm scheduler; has an RTX 3080 GPU; /ngs is mounted so the conda
# env and input files are already present).
#   STEP 1: imputation (probability space, self-verified)  -> output_imp_test/
#   STEP 2: recontextualization (logit space, diagnostic)  -> output_recon_test/
#
# Every line of stdout+stderr is tee'd to a timestamped log:
#     /ngs/iflores/andrea/logs/shap_tests_<timestamp>.log
# so a crash is always traceable.
#
# Run it directly (NOT via a scheduler). To survive SSH disconnects use tmux:
#     tmux new -s shap
#     bash /ngs/iflores/andrea/run_shap_tests.sh
#     # Ctrl-b then d  to detach;  tmux attach -t shap  to come back

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
# carmela has conda available as a shell function and the tabpfn2 env on /ngs.
source /ngs/software/conda/etc/profile.d/conda.sh 2>/dev/null || true
conda activate /ngs/software/conda/envs/tabpfn2 \
  || micromamba activate /ngs/software/conda/envs/tabpfn2 \
  || { echo "ERROR: could not activate tabpfn2 env"; exit 1; }
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
  --n-jobs      4  \
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
  --n-jobs      4  \
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
