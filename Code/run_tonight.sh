#!/usr/bin/env bash
# =============================================================================
# run_tonight.sh -- one-shot finalisation batch (launch once, walk away).
#
# Launch from the project root, wrapped in caffeinate so the Mac stays awake:
#
#   caffeinate -dims bash Code/run_tonight.sh 2>&1 | tee -a Results/tonight_launch.log
#
# Arguments:
#   $1  worker count            (default 8)
#   $2  resume-from stage id    (default 0 = run everything). E.g. "2" skips
#       the unit gates and the dgp2 pilot gate after a verified earlier pass.
#
# Stages (cheapest / most tex-blocking first; ids used by the resume argument):
#   0  unit gates (fatal on failure)
#   1  PILOT GATE for the second DGP (fatal on failure; ~5-30 min)
#   2  MD&A restart dispersion at K = 50 (three single-start refits, ~45 min)
#   3  E4 correct-prior size arm (alpha = 0.5 = generating value; size only)
#   4  Study-I selection on 100 replicates (evaluation-only, cached fits)
#   5  second DGP, Study-I scope: E1 selection + E2 coverage (new fits, W=10k)
#   6  MD&A sensitivity: c = 0.5, c = 2, truncated grid (evaluation-only)
#   7  post-processing + figure/table regeneration
#
# Stages 2-7 are RESILIENT: a failure is logged and the chain continues, so one
# broken arm cannot take the rest of the night down. The script exits non-zero
# if anything failed and prints the failed list at the end.
#
# Rough ETA on ~8 workers: 8-11 h nominal; the shipped run of 1-2 Sep 2026
# took 16 h (selreps alone ~8 h). Every stage appends to
# Results/tonight_<stamp>.log and prints its own duration.
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")/.."

W="${1:-8}"
FROM="${2:-0}"
STAMP=$(date +%Y%m%d_%H%M)
LOG="Results/tonight_${STAMP}.log"
mkdir -p Results
FAILED=()

# --- paper-design override blocks -------------------------------------------
BASE_DGP="K_true=40 J_train=1000 W=5000 K_grid=10:100:10 fit_method=WarpLDA beta=0.01"
DGP2="K_true=20 J_train=1000 W=10000 L=500 K_grid=5:50:5 fit_method=WarpLDA alpha=0.1 beta=0.01"

say() { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOG"; }

# fatal stage: abort the whole chain on failure (gates / pilot only)
stage_fatal() { # stage_fatal <id> <name> <cmd...>
  local id=$1 name=$2; shift 2
  if [ "$FROM" -gt "$id" ]; then say "---------- SKIP  [$id] $name (resume from $FROM)"; return 0; fi
  local t0=$SECONDS
  say "========== START [$id] $name =========="
  if ! "$@" >>"$LOG" 2>&1; then
    say "========== FATAL [$id] $name failed -- aborting chain =========="
    exit 1
  fi
  say "========== DONE  [$id] $name in $(( (SECONDS - t0) / 60 )) min =========="
}

# resilient stage: log the failure, keep going
stage() { # stage <id> <name> <cmd...>
  local id=$1 name=$2; shift 2
  if [ "$FROM" -gt "$id" ]; then say "---------- SKIP  [$id] $name (resume from $FROM)"; return 0; fi
  local t0=$SECONDS
  say "========== START [$id] $name =========="
  if "$@" >>"$LOG" 2>&1; then
    say "========== DONE  [$id] $name in $(( (SECONDS - t0) / 60 )) min =========="
  else
    say "========== FAIL  [$id] $name after $(( (SECONDS - t0) / 60 )) min (chain continues) =========="
    FAILED+=("[$id] $name")
  fi
}

# --- preflight ---------------------------------------------------------------
say "run_tonight.sh | workers=$W | resume-from=$FROM | log=$LOG"
say "disk: $(df -h . | tail -1 | awk '{print $4 " free"}')"
[ -d Data/FITS ] || { say "FATAL: Data/FITS missing"; exit 1; }
[ -f Data/MDNA/mdna_prep_2015_2016.qs2 ] || { say "FATAL: MD&A prep missing"; exit 1; }
say "ETA: gates ~5m | pilot ~30m | restarts ~45m | E4 prior ~1-2h | selreps ~1-2h | dgp2 ~2-4h | MD&A 3x ~5h"

# --- 0. unit gates -----------------------------------------------------------
stage_fatal 0 "unit gates (smoke)" Rscript Code/tests_unit.R smoke

# --- 1. pilot gate for the second DGP ---------------------------------------
stage_fatal 1 "dgp2 PILOT gate" Rscript Code/run_E1.R pilot "$W" $DGP2 \
      S=2 J_eval=1000 n_starts=3 label=dgp2PILOT

# --- 2. MD&A restart dispersion (cheap; unblocks a main-text TODO) -----------
stage 2 "MD&A restart dispersion (K=50)" Rscript Code/run_mdna_restarts.R \
      workers="$W" K=50

# --- 3. E4 correct-prior size arm -------------------------------------------
stage 3 "E4 correct-prior size arm" Rscript Code/run_E4.R full "$W" $BASE_DGP \
      alpha=0.5 S_train=10 R_null=500 R_power=0 alternatives=none \
      J_ev=250 J_center=10000 K_test=10,40 n_starts=3 label=E4prior05

# --- 4. Study-I selection on 100 replicates (evaluation-only) ----------------
stage 4 "E1 selreps (10x10, Jev=5000)" Rscript Code/run_E1_selreps.R full "$W" \
      $BASE_DGP alpha=0.1 J_eval=5000 S=10 n_starts=3 R_sel=10 label=selreps

# --- 5. second DGP: Study-I scope -------------------------------------------
stage 5 "dgp2 E1 (selection)" Rscript Code/run_E1.R full "$W" $DGP2 \
      J_eval=5000 S=10 n_starts=3 label=dgp2
stage 5 "dgp2 E2 (coverage)" Rscript Code/run_E2.R full "$W" $DGP2 \
      S_train=10 R_eval=10 J_ev_grid=100,250,500 J_truth=2000 \
      K_test=10,20,30 n_starts=3 label=dgp2

# --- 6. MD&A sensitivity (evaluation-only, cached fits) ----------------------
stage 6 "MD&A c=0.5" Rscript Code/run_mdna.R workers="$W" K_grid=10:200:10 \
      tests_K=hat c_part=0.5 out_suffix=_c05
stage 6 "MD&A c=2"   Rscript Code/run_mdna.R workers="$W" K_grid=10:200:10 \
      tests_K=hat c_part=2 out_suffix=_c2
stage 6 "MD&A grid 10:100" Rscript Code/run_mdna.R workers="$W" \
      K_grid=10:100:10 tests_K=hat out_suffix=_g100

# --- 7. post-processing + regeneration ---------------------------------------
stage 7 "postprocess E4 (paper object)" Rscript Code/postprocess_e4.R \
      file=Data/E4/e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2
stage 7 "postprocess E4 (prior arm, newest object)" Rscript Code/postprocess_e4.R
stage 7 "postprocess E2 gap" Rscript Code/postprocess_e2_gap.R
stage 7 "figures (paper label)" Rscript Code/make_figures.R full $BASE_DGP \
      alpha=0.1 label=warptestFULL_NULL
stage 7 "tables (paper label)" Rscript Code/make_tables.R full $BASE_DGP \
      alpha=0.1 label=warptestFULL_NULL
stage 7 "figures (dgp2; best-effort)" Rscript Code/make_figures.R full $DGP2 \
      exp=E1,E2 label=dgp2

# --- summary -----------------------------------------------------------------
if [ "${#FAILED[@]}" -gt 0 ]; then
  say "=== CHAIN FINISHED WITH FAILURES: ${FAILED[*]} -- see $LOG ==="
  exit 1
fi
say "=== ALL STAGES DONE — results current; summaries in $LOG ==="
say "next: review the log, then the tex fold-in (Phase C)."
