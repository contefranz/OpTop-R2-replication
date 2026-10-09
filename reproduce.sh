#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# reproduce.sh -- entry points for the replication lanes (see README.md).
#
#   ./reproduce.sh test                unit gates U1-U25 (~10 s; run first)
#   ./reproduce.sh fetch WHAT [--local DIR]
#                                      download (or take from DIR) the deposited data,
#                                      verify every checksum and unpack it in place:
#                                      WHAT = results (~270 MB) | fits (2.7 GB) |
#                                             corpus (178 MB) | all
#   ./reproduce.sh paper               DEFAULT. Every table and figure of the paper and
#                                      the Supplementary Material, from the shipped caches
#                                      and the deposited result objects (~2 min; needs
#                                      `fetch results`). Overwrites Results/ in place.
#   ./reproduce.sh stage3 [N]          the computations that read the fit cache: MD&A
#                                      completion-residual tests (~4 min) and the per-fit
#                                      completion reference curves (~40 min on N=8
#                                      workers); needs `fetch results` and `fetch fits`
#   ./reproduce.sh revision [N]        the evaluation-only revision batch behind the *_rev2
#                                      result objects (~28 h on N=8 workers; needs `fetch
#                                      fits`); completed stages whose inputs are unchanged
#                                      are verified and skipped (README.md, Section 5)
#   ./reproduce.sh legacy              the pre-revision exhibit trees (2 Sep 2026 Lane A)
#   ./reproduce.sh smoke               tiny end-to-end run of the simulation pipeline
#   ./reproduce.sh full-sims [N]       refit every simulation design (days; refits are
#                                      statistically equivalent, not bit-identical)
#   ./reproduce.sh full-mdna [N]       MD&A corpus prep + estimation (needs `fetch corpus`)
#   ./reproduce.sh robustness [N]      the September 2026 robustness arms (refits)
# ------------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
# a project-local library created by `Rscript install.R lib=.Rlib` takes precedence
if [ -d .Rlib ]; then export R_LIBS_USER="$PWD/.Rlib"; fi

# Zenodo data record holding the large objects (README.md, Section 4)
ZENODO_RECORD="${ZENODO_RECORD:-RECORD_ID_TO_BE_FILLED}"
ZENODO_URL="https://zenodo.org/records/${ZENODO_RECORD}/files"

D=(K_true=40 W=5000 K_grid=10:100:10 fit_method=WarpLDA alpha=0.1 beta=0.01)
SIM_REPS=(J_eval=5000 S=10 R_eval=10 J_truth=2000 J_ev=250 R_power=250 J_center=10000 K_test=10,40)
DGP2=(K_true=20 J_train=1000 W=10000 L=500 K_grid=5:50:5 fit_method=WarpLDA alpha=0.1 beta=0.01)
E4_REV=Data/E4/e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2
E2_REV40=Data/E2/e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2
E2_REV20=Data/E2/e2_results_E2_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01_rev2.qs2

run() { echo "+ $*"; "$@"; }

need() {   # need <manifest> <fetch target>: every file listed in the tracked manifest must exist
  local missing
  missing=$(awk '{ $1 = ""; sub(/^ +/, ""); print }' "$1" | while IFS= read -r f; do [ -f "$f" ] || echo "$f"; done | wc -l | tr -d ' ')
  if [ "$missing" != "0" ]; then
    echo "$missing deposited file(s) listed in $1 are missing. Run first:  ./reproduce.sh fetch $2" >&2
    exit 1
  fi
}
need_results() { need Data/MANIFEST-result_objects.sha256 results; }
need_fits()    { need Data/MANIFEST-fit_cache.sha256 fits; }

fetch_one() {   # fetch_one <archive> [local dir]: checksum of the archive, then of every member
  local name="$1" src="${2:-}" tmp manifest
  manifest="Data/MANIFEST-${name%.tar}.sha256"     # tracked in git: the expected contents
  tmp="$(mktemp -d)"
  if [ -n "$src" ]; then
    cp "$src/$name" "$src/SHA256SUMS" "$tmp/"
  else
    if [ "$ZENODO_RECORD" = "RECORD_ID_TO_BE_FILLED" ]; then
      echo "The Zenodo record id is not set yet (ZENODO_RECORD in reproduce.sh)." >&2; exit 1
    fi
    curl -fL --retry 3 -o "$tmp/SHA256SUMS" "$ZENODO_URL/SHA256SUMS?download=1"
    curl -fL --retry 3 -o "$tmp/$name" "$ZENODO_URL/$name?download=1"
  fi
  (cd "$tmp" && grep " $name\$" SHA256SUMS | shasum -a 256 -c -)
  tar -xf "$tmp/$name"
  shasum -a 256 -c --quiet "$manifest" && echo "verified: $name ($(wc -l < "$manifest" | tr -d ' ') files)"
  rm -rf "$tmp"
}

case "${1:-paper}" in
  test)
    run Rscript Code/tests_unit.R ;;

  fetch)
    what="${2:-}"; [ -n "$what" ] || { echo "usage: ./reproduce.sh fetch {results|fits|corpus|all} [--local DIR]" >&2; exit 1; }
    src=""; if [ "${3:-}" = "--local" ]; then src="${4:?--local needs a directory}"; fi
    case "$what" in
      results) fetch_one result_objects.tar "$src" ;;
      fits)    fetch_one fit_cache.tar "$src" ;;
      corpus)  fetch_one mdna_corpus.tar "$src" ;;
      all)     for a in result_objects.tar fit_cache.tar mdna_corpus.tar; do fetch_one "$a" "$src"; done ;;
      *) echo "unknown fetch target: $what" >&2; exit 1 ;;
    esac ;;

  paper)
    need_results
    # revision post-processing: manuscript-facing csv and LaTeX bodies (Results/csv, Results/tex)
    run Rscript Code/postprocess_rev2.R
    run Rscript Code/postprocess_revision.R
    run Rscript Code/make_revision_tex.R
    run Rscript Code/make_revision_figures.R
    run Rscript Code/postprocess_e4.R     file=$E4_REV
    run Rscript Code/postprocess_e2_gap.R file=$E2_REV40
    run Rscript Code/postprocess_e2_gap.R file=$E2_REV20
    # supplement figures (S1-S11) and the pre-revision tables still cited (S2, S8, S9, S14, S15)
    run Rscript Code/make_figures.R full exp=E1,E2,E4,E6 "${D[@]}" in_suffix=_rev2 label=warptestFULL_NULL_rev2
    run Rscript Code/make_figures.R full exp=E3,E5 "${D[@]}" label=warptestFULL_NULL
    run Rscript Code/make_tables.R  full "${D[@]}" label=warptestFULL_NULL
    # empirical application (Figure S12 and the MD&A tables)
    run Rscript Code/make_mdna_outputs.R y1=2015 y2=2016 in_suffix=_rev2
    # main Figures 1-4 and the Stage 3 computations that need no fit cache (Results/manuscript)
    run Rscript Code/manuscript/make_main_figures.R
    run Rscript Code/manuscript/stage3_e2_falsecert.R
    run Rscript Code/manuscript/stage3_e1_selreps.R
    run Rscript Code/manuscript/stage3_mdna_subsamples.R
    run Rscript Code/manuscript/stage3_mdna_H_bootstrap.R
    run Rscript Code/manuscript/s64_wald_decomposition.R
    run Rscript Code/manuscript/stage3_completion_paired_se.R
    run Rscript Code/manuscript/stage3_make_tables.R
    echo "Done. Manuscript figure map: Results/manuscript_files.csv; exhibit map: README.md, Section 6." ;;

  stage3)
    need_results; need_fits
    run env OPTOP_NO_FIT=1 Rscript Code/manuscript/stage3_mdna_completion_moments.R
    run env OPTOP_NO_FIT=1 Rscript Code/manuscript/stage3_completion_reference.R full "${2:-8}"
    run Rscript Code/manuscript/stage3_completion_paired_se.R
    run Rscript Code/manuscript/stage3_make_tables.R ;;

  revision)
    need_fits
    run bash Code/run_revision_batch.sh "${2:-8}" A _rev2
    run bash Code/run_revision_batch.sh "${2:-8}" B _rev2
    run bash Code/run_revision_batch.sh "${2:-8}" final _rev2
    echo "Batch complete; now run  ./reproduce.sh paper" ;;

  legacy)
    run Rscript Code/make_figures.R full "${D[@]}" label=warptestFULL_NULL
    run Rscript Code/make_tables.R  full "${D[@]}" label=warptestFULL_NULL
    run Rscript Code/postprocess_e4.R     file=Data/E4/e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2
    run Rscript Code/postprocess_e2_gap.R file=Data/E2/e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2
    run Rscript Code/postprocess_e2_gap.R file=Data/E2/e2_results_E2_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01.qs2
    run Rscript Code/make_figures.R full "${DGP2[@]}" exp=E1,E2 label=dgp2
    run Rscript Code/make_mdna_outputs.R y1=2015 y2=2016 ;;

  smoke)      run Rscript Code/make_all.R smoke ;;
  full-sims)  run Rscript Code/make_all.R full "${2:-11}" "${D[@]}" "${SIM_REPS[@]}" label=warptestFULL_NULL ;;
  full-mdna)
    [ -f Data/corpus_item7/item7_2016_2017.qs2 ] || { echo "Run first:  ./reproduce.sh fetch corpus" >&2; exit 1; }
    run Rscript Code/prep_mdna.R y1=2015 y2=2016 workers="${2:-6}"
    run Rscript Code/run_mdna.R  y1=2015 y2=2016 workers="${2:-6}" K_grid=10:200:10 refine_span=0
    run Rscript Code/make_mdna_outputs.R y1=2015 y2=2016 ;;
  robustness) run bash Code/run_tonight.sh "${2:-8}" ;;
  *) echo "usage: ./reproduce.sh {test|fetch WHAT|paper|stage3 [N]|revision [N]|legacy|smoke|full-sims [N]|full-mdna [N]|robustness [N]}" >&2; exit 1 ;;
esac
