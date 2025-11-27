#!/usr/bin/env bash
# Run compare_JEC_SFs.py (must succeed), then run two independent tasks in parallel:
#  - python3 check_L1JetSFs.py
#  - cd ../makeLUTs && python3 updateSFPtEtaBins.py
#
# All stdout/stderr is written to out.txt (and also echoed to the terminal).
# Exit only if compare_JEC_SFs.py fails. If either parallel job fails, report which one failed
# but do NOT exit non-zero because of those failures.
#
# Usage:
#   chmod +x afterBDT.sh
#   ./afterBDT.sh
set -u

LOG="out.txt"
# Start fresh log
: > "$LOG"
# Tee all stdout/stderr to the log while still printing to the terminal
exec > >(tee -a "$LOG") 2>&1

timestamp() { date +%FT%T; }

echo "[$(timestamp)] Starting afterBDT.sh"

echo "[$(timestamp)] Running: python3 compare_JEC_SFs.py"
if python3 compare_JEC_SFs.py; then
  echo "[$(timestamp)] compare_JEC_SFs.py completed SUCCESSFULLY"
else
  rc=$?
  echo "[$(timestamp)] ERROR: compare_JEC_SFs.py failed with exit code $rc" >&2
  echo "[$(timestamp)] Exiting because compare_JEC_SFs.py failed"
  exit $rc
fi

echo "[$(timestamp)] Launching two jobs in parallel..."

# Launch job A: check_L1JetSFs.py in current directory
(
  echo "[$(timestamp)] [check_L1JetSFs.py] Starting"
  python3 check_L1JetSFs.py
  rc=$?
  echo "[$(timestamp)] [check_L1JetSFs.py] Exited with code $rc"
  exit $rc
) &
PID_A=$!
echo "[$(timestamp)] Started check_L1JetSFs.py (PID $PID_A)"

# Launch job B: updateSFPtEtaBins.py inside ../makeLUTs
(
  echo "[$(timestamp)] [../makeLUTs/updateSFPtEtaBins.py] Starting"
  cd ../makeLUTs || { echo "[$(timestamp)] ERROR: cannot cd to ../makeLUTs"; exit 127; }
  python3 updateSFPtEtaBins.py
  rc=$?
  echo "[$(timestamp)] [../makeLUTs/updateSFPtEtaBins.py] Exited with code $rc"
  exit $rc
) &
PID_B=$!
echo "[$(timestamp)] Started updateSFPtEtaBins.py in ../makeLUTs (PID $PID_B)"

# Wait for both and report results. Do NOT exit non-zero based on these results.
ANY_FAILED=0

wait "$PID_A"
RC_A=$?
if [ $RC_A -eq 0 ]; then
  echo "[$(timestamp)] check_L1JetSFs.py (PID $PID_A) finished SUCCESS"
else
  echo "[$(timestamp)] ERROR: check_L1JetSFs.py (PID $PID_A) failed with exit code $RC_A"
  ANY_FAILED=1
fi

wait "$PID_B"
RC_B=$?
if [ $RC_B -eq 0 ]; then
  echo "[$(timestamp)] updateSFPtEtaBins.py (PID $PID_B) finished SUCCESS"
else
  echo "[$(timestamp)] ERROR: updateSFPtEtaBins.py (PID $PID_B) failed with exit code $RC_B"
  ANY_FAILED=1
fi

if [ $ANY_FAILED -ne 0 ]; then
  echo "[$(timestamp)] One or more parallel jobs failed. See $LOG for full output."
else
  echo "[$(timestamp)] All parallel jobs finished successfully."
fi

echo "[$(timestamp)] afterBDT.sh finished (exiting 0 because compare_JEC_SFs.py succeeded)"
exit 0