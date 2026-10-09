#!/usr/bin/env bash
# =============================================================================
# revision_smoke.sh -- every driver of the revision batch, end to end, on toy
# corpora, in a few minutes. The ONLY batch stage allowed to estimate models
# (toy ones, as the unit suite does). It exercises, in this order:
#   1. each simulation driver with out_suffix (checkpoints written);
#   2. the RESUME path: E1 again under the same suffix must reuse every
#      checkpoint and reproduce the object bit for bit;
#   3. cfg_from= + OPTOP_NO_FIT=1: re-scoring from a stored configuration finds
#      every fit in the cache;
#   4. the MD&A driver on the 300-document sample with c != delta and an explicit
#      tests_K list, the restart script, and both diagnostics modes.
# Usage: bash Code/revision_smoke.sh [workers=2] [suffix=_rev2]
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."
W="${1:-2}"; SFX="${2:-_rev2}_smoke"
case "$SFX" in *_smoke) ;; *) echo "refusing: smoke suffix must end in _smoke"; exit 2;; esac

# a smoke run always starts from scratch: checkpoints are tied to the code that
# wrote them, so stale ones from an earlier version of the code must go
rm -rf "Data/Checkpoints/${SFX}" "Data/Checkpoints/${SFX}_replay"
say() { printf '\n[smoke %s] %s\n' "$(date +%H:%M:%S)" "$*"; }

for DRIVER in E1 E1_selreps E2 E4 E6; do
  say "run_${DRIVER}.R"
  extra=(); [ "$DRIVER" = "E1_selreps" ] && extra=("R_sel=2")
  # ${arr[@]+...}: an empty array is "unbound" under set -u in macOS bash 3.2
  Rscript "Code/run_${DRIVER}.R" smoke "$W" ${extra[@]+"${extra[@]}"} "out_suffix=$SFX"
done

E1_OBJ=$(Rscript --vanilla -e 'suppressMessages(source("Code/R/source_all.R")); source("Code/config/configs.R"); cat(file.path("Data/E1", paste0("e1_results_", run_tag(get_config("E1", "smoke")), commandArgs(TRUE)[1], ".qs2")))' "$SFX" | tail -1)
[ -f "$E1_OBJ" ] || { echo "smoke: E1 object not found: $E1_OBJ"; exit 1; }

say "resume path: E1 again under the same suffix (checkpoints must be reused)"
H0=$(shasum -a 256 "$E1_OBJ" | cut -d' ' -f1); cp "$E1_OBJ" "${E1_OBJ}.first"
Rscript Code/run_E1.R smoke "$W" "out_suffix=$SFX"
Rscript --vanilla -e 'suppressMessages(library(qs2)); a <- qs_read(commandArgs(TRUE)[1]); b <- qs_read(commandArgs(TRUE)[2]); attr(a, "run_meta") <- attr(b, "run_meta") <- NULL; if (!isTRUE(all.equal(a, b, tolerance = 0))) stop("resumed E1 object differs from the first run"); cat("resume: identical object\n")' "${E1_OBJ}.first" "$E1_OBJ"
rm -f "${E1_OBJ}.first"

say "cfg_from + no-fit guard: re-score the stored configuration (cache hits only)"
OPTOP_NO_FIT=1 Rscript Code/run_E1.R smoke "$W" "cfg_from=$E1_OBJ" "out_suffix=${SFX}_replay"

say "run_mdna.R on the 300-document sample (c = 0.5, delta = 1, explicit tests_K)"
Rscript Code/run_mdna.R sample_n=300 "workers=$W" K_grid=10:30:10 refine_span=0 \
        tests_K=10,30 c_part=0.5 min_null=1 "out_suffix=$SFX"

say "run_mdna_restarts.R on the sample (support common to grid and restarts)"
Rscript Code/run_mdna_restarts.R sample_n=300 "workers=$W" K=20 common=2 \
        K_grid=10:30:10 "out_suffix=$SFX"

say "run_revision_diagnostics.R (toy sources)"
Rscript Code/run_revision_diagnostics.R "$W" mode=unbinned "sources=$E1_OBJ" \
        mdna_sample_n=300 mdna_grid=10:30:10 "out_suffix=$SFX"
E5_OBJ=$(Rscript --vanilla -e 'suppressMessages(source("Code/R/source_all.R")); source("Code/config/configs.R"); cat(file.path("Data/E5", paste0("e5_results_", run_tag(get_config("E5", "smoke")), ".qs2")))' | tail -1)
if [ ! -f "$E5_OBJ" ]; then say "run_E5.R (smoke object absent)"; Rscript Code/run_E5.R smoke "$W"; fi
Rscript Code/run_revision_diagnostics.R "$W" mode=e5 "input=$E5_OBJ" "out_suffix=$SFX"

say "smoke complete"
