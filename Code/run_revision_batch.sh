#!/usr/bin/env bash
# =============================================================================
# run_revision_batch.sh -- launcher of the evaluation-only revision batch
# (Code/revision_batch.py holds the stage table and the completion logic).
#
#   bash Code/run_revision_batch.sh [workers=8] [select=all] [suffix=_rev2] [--dry-run]
#
#   select   A      night 1 : MD&A (2), MD&A designs (3), E2 (4), E1 (5), E6 (6),
#                             E5 resolution (9), unbinned comparator (10), restarts (11)
#            B      nights 2-3: E4 both prior arms (7), E1 10 x 10 selection (8)
#            final  post-processing (12) and acceptance checks (13)
#            all | an id range "2-6" | ids or keys "mdna,e6,7"
#
# Typical use (each line resumes verified work and can be re-issued at will):
#   bash Code/run_revision_batch.sh 8 A     2>&1 | tee -a Results/revision_rev2/launch.log
#   bash Code/run_revision_batch.sh 8 B     2>&1 | tee -a Results/revision_rev2/launch.log
#   bash Code/run_revision_batch.sh 8 final 2>&1 | tee -a Results/revision_rev2/launch.log
#
# Unit gates and the toy smoke run are prerequisites: they are (re)run
# automatically whenever the scoring code or the package changed. Stage status:
# Results/revision<suffix>/stages.json, summary.json and one log per stage.
# The machine is kept awake with caffeinate when available.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
W="${1:-8}"; SEL="${2:-all}"; SFX="${3:-_rev2}"
shift $(( $# < 3 ? $# : 3 ))
mkdir -p "Results/revision${SFX}"
if command -v caffeinate >/dev/null 2>&1; then
  exec caffeinate -dims python3 Code/revision_batch.py --workers "$W" --select "$SEL" --suffix "$SFX" "$@"
else
  exec python3 Code/revision_batch.py --workers "$W" --select "$SEL" --suffix "$SFX" "$@"
fi
