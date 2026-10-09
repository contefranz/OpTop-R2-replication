# Code/manuscript — main-text figures and the Stage 3 computations

Scripts written for the final manuscript (24 September - 9 October 2026). They read saved
results only, never fit a model, and write to `Results/manuscript/` (figures and sidecar
csv) and `Results/manuscript/stage3/` (Stage 3 outputs). Run them from the package root;
`./reproduce.sh paper` runs those that need no fit cache, `./reproduce.sh stage3` the others.

| Script | Produces | Manuscript use | Needs | Time |
|---|---|---|---|---|
| `make_main_figures.R` | `F1_selection_support.pdf`, `gap_decomposition.pdf`, `F2_mdna_fit.pdf`, `residual_mass.pdf` + 13 sidecar csv | Figures 1-4; Table 2 "Macro maximised" rows | E1 and MD&A `_rev2` objects | 5 s |
| `stage3_e2_falsecert.R` | `stage3_e2_falsecert.csv`, `stage3_e2_reference_curves.csv` | Table S5, Section 5.1 | E2 `_rev2` object | 1 s |
| `stage3_e1_selreps.R` | `stage3_e1_falsecert.csv`, `stage3_e1_table2_mcse.csv` | Table S5 (J_ev = 5,000), Table 2 Monte Carlo SEs | E1 selection replicates; previous script | 1 s |
| `stage3_mdna_subsamples.R` | `stage3_mdna_subsample_{selections,levels}.csv`, `stage3_mdna_gap_cluster.csv` | Table S28 present-firm rows; clustered gap intervals | MD&A `_rev2` object, prep | 7 s |
| `stage3_mdna_H_bootstrap.R` | `stage3_mdna_H_bootstrap.csv` | bootstrap interval for H (Section 6.3, S6.4) | `residual_cluster_sums.csv` (from `make_main_figures.R`) | 3 s |
| `s64_wald_decomposition.R` | `stage3_mdna_s64_wald_subtests.csv` | S6.4 sub-tests; Section 6.3 t = 3.2, p = 0.08 | sidecar csv | < 1 s |
| `stage3_completion_reference.R {smoke,pilot,full} [workers]` | `stage3_completion_reference_curves[_mode].csv`, `stage3_completion_falsecert.csv`, `Data/E2/e2_completion_truth_rev2[_mode].qs2` | Table S5 completion block, "880 certificates" | **fit cache** | 16 s / 5 min / 40 min (8 workers) |
| `stage3_completion_paired_se.R` | `stage3_completion_paired_se.csv` | S5.1 and Table S5 notes (fit 7: 0.00086) | completion truth object | 2 s |
| `stage3_mdna_completion_moments.R` | `stage3_mdna_completion_{tests,strata,masses}.csv`, `Data/MDNA/stage3_mdna_completion_G.qs2` | Table S25; Sections 4.1, 6.3 | **fit cache** (MD&A K = 50, 180) | 4 min (8 workers) |
| `stage3_make_tables.R` | `stage3_falsecert_{fits,table}.csv`, `tabS5_rows.tex`, `tabS25_rows.tex`, `present_rows.tex` | LaTeX rows of Tables S5, S25, S28 | the Stage 3 csv | < 1 s |

Order: `stage3_e2_falsecert.R` before `stage3_e1_selreps.R`; `make_main_figures.R` before
`stage3_mdna_H_bootstrap.R`; the completion scripts before `stage3_make_tables.R`.
Table S5 as printed merges the reconstruction rows of `tabS5_rows.tex` with the completion
columns and block from `stage3_completion_reference_curves.csv` and
`stage3_completion_falsecert.csv`.

Every script except `s64_wald_decomposition.R` and `stage3_completion_paired_se.R` is the
script the authors ran on 7-8 October 2026, with input and output paths adapted to the
package; on 9 October 2026 each one, run in a clean copy of the package, reproduced its
saved outputs byte for byte (smoke and pilot modes for the completion reference curves).
