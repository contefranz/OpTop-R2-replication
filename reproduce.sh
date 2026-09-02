#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# reproduce.sh -- entry points for the replication lanes (see README.md).
#
#   ./reproduce.sh test               unit gates (run first)
#   ./reproduce.sh exhibits           Lane A: all figures + tables from the
#                                     shipped caches, incl. post-processing
#                                     and the second-DGP figures (~minutes,
#                                     no refitting)
#   ./reproduce.sh smoke              tiny end-to-end run of the simulation
#                                     pipeline (~2 min; exercises all code paths)
#   ./reproduce.sh full-sims [N]      Lane B: definitive simulation run on N
#                                     workers (default 11; ~1 day wall-clock)
#   ./reproduce.sh full-mdna [N]      Lane B: MD&A prep + estimation + outputs
#                                     on N workers (default 6; ~7 h wall-clock)
#   ./reproduce.sh robustness [N]     Lane C: the September 2026 robustness
#                                     arms (size-adjusted power, correct-prior
#                                     size, replicated selection, second DGP,
#                                     MD&A c/grid sensitivity, restart
#                                     dispersion); ~16 h on 8 workers
# ------------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"

SIM_DESIGN="K_true=40 W=5000 K_grid=10:100:10 fit_method=WarpLDA alpha=0.1 beta=0.01"
SIM_REPS="J_eval=5000 S=10 R_eval=10 J_truth=2000 J_ev=250 R_power=250 J_center=10000 K_test=10,40"
DGP2_DESIGN="K_true=20 J_train=1000 W=10000 L=500 K_grid=5:50:5 fit_method=WarpLDA alpha=0.1 beta=0.01"
E4_OBJ="e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"
E2_OBJ="e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"
E2_OBJ_DGP2="e2_results_E2_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01.qs2"

case "${1:-exhibits}" in
  test)
    Rscript Code/tests_unit.R
    ;;
  exhibits)
    Rscript Code/make_figures.R full $SIM_DESIGN label=warptestFULL_NULL
    Rscript Code/make_tables.R  full $SIM_DESIGN label=warptestFULL_NULL
    Rscript Code/postprocess_e4.R     file=Data/E4/$E4_OBJ
    Rscript Code/postprocess_e2_gap.R file=Data/E2/$E2_OBJ
    Rscript Code/postprocess_e2_gap.R file=Data/E2/$E2_OBJ_DGP2
    Rscript Code/make_figures.R full $DGP2_DESIGN exp=E1,E2 label=dgp2
    Rscript Code/make_mdna_outputs.R y1=2015 y2=2016
    echo "Done. Figures: Results/Figures/{warptestFULL_NULL,dgp2,mdna_2015_2016}/"
    echo "      Tables : Results/{csv,tex,xlsx}/"
    ;;
  smoke)
    Rscript Code/make_all.R smoke
    ;;
  full-sims)
    Rscript Code/make_all.R full "${2:-11}" $SIM_DESIGN $SIM_REPS label=warptestFULL_NULL
    ;;
  full-mdna)
    Rscript Code/prep_mdna.R y1=2015 y2=2016 workers="${2:-6}"
    Rscript Code/run_mdna.R  y1=2015 y2=2016 workers="${2:-6}" K_grid=10:200:10 refine_span=0
    Rscript Code/make_mdna_outputs.R y1=2015 y2=2016
    ;;
  robustness)
    bash Code/run_tonight.sh "${2:-8}"
    ;;
  *)
    echo "usage: ./reproduce.sh {test|exhibits|smoke|full-sims [N]|full-mdna [N]|robustness [N]}" >&2
    exit 1
    ;;
esac
