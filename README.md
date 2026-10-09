# Replication Package

**Goodness-of-Fit Indices and Diagnostics for Topic Models**
Craig M. Lewis (Vanderbilt University) and Francesco Grossetti (Bocconi University)
Submitted to the *Journal of the Royal Statistical Society, Series B*. Contact: francesco.grossetti@unibocconi.it

- Code: https://github.com/contefranz/OpTop-R2-replication (archived on Zenodo: DOI *to be added at release*)
- Data objects: Zenodo data record, DOI *to be added at release* (Section 3)
- Licences: code MIT (`LICENSE`); data and results CC0-1.0 (`LICENSE-DATA`)

This package reproduces **every table, figure and reported number** of the paper (Tables 1-4, Figures 1-4) and of its Supplementary Material (Tables S1-S29, Figures S1-S13): the four simulation studies, the MD&A empirical application, the evaluation-only revision analyses and the Stage 3 computations. The harmonised support and the discrepancy indices are computed with the R package [`OpTop`](https://github.com/contefranz/OpTop) 0.20.1; the held-out protocols, the certification rules (total-gain and adjacent), the firm-clustered inference, the residual diagnostics and the held-out word-level baseline are implemented here, in `Code/R/`.

All exhibits are regenerated **exactly** from cached results: on 9 October 2026 a clean copy of this package regenerated every csv and LaTeX file byte for byte and every figure identically up to PDF creation dates (Section 8). Nothing is re-estimated in the default lane.

---

## Quick start

```sh
Rscript install.R lib=.Rlib     # exact environment in a project-local library (Section 2)
./reproduce.sh test             # unit gates U1-U25 (~10 s)
./reproduce.sh fetch results    # result objects from the Zenodo data record (259 MB, checksummed)
./reproduce.sh paper            # every exhibit of the paper and the supplement (~1.5 min)
```

`./reproduce.sh paper` overwrites `Results/` in place; `git status` afterwards shows only PDFs (new creation dates) and the two Excel workbooks (new timestamps). `Results/manuscript_files.csv` maps every figure file of the manuscript (`Figures/<name>.pdf`) to its file in this package, with its sha256.

## 1. Contents

```
README.md, LICENSE, LICENSE-DATA, CITATION.cff, .zenodo.json
install.R              exact environment: R 4.6.1, CRAN snapshot 2026-07-01, OpTop/NLPstudio commits
reproduce.sh           entry points for all lanes (Section 4)
sessionInfo.txt        the environment behind every shipped result, timeline, per-object stamps
.here                  project-root anchor (paths resolve against it; do not delete)
Code/
  R/                   analysis modules: I/O, DGP, fitting, held-out scoring, inference
                       (certificates, clustered variances), moment tests, comparators, revision
  config/configs.R     design parameters, profiles (smoke/pilot/full), seeds
  run_E1.R..run_E6.R   simulation drivers (cache-backed, resumable); run_E1_selreps.R
  prep_mdna.R, run_mdna.R, run_mdna_restarts.R      empirical application
  make_figures.R, make_tables.R, make_mdna_outputs.R exhibit builders (read saved results only)
  postprocess_e4.R, postprocess_e2_gap.R            post-processing of saved E4/E2 objects
  rescore_revision.R, postprocess_revision.R, make_revision_tex.R, make_revision_figures.R
                       September 2026 revision, phase rev1 (MD&A rescoring and its exhibits)
  revision_batch.py, run_revision_batch.sh, revision_smoke.sh, run_revision_diagnostics.R,
  finalize_revision.R, postprocess_rev2.R
                       evaluation-only revision batch (phase rev2), its acceptance and summaries
  manuscript/          main Figures 1-4 and the Stage 3 computations (Code/manuscript/README.md)
  make_all.R, run_tonight.sh, list_runs.R, tests_unit.R
  README.md            development manual of the pipeline (CLI grammar, engines, design notes)
Data/
  E1/..E6/, MDNA/      pre-revision result caches (tracked); deposited objects unpack here
  MANIFEST-*.sha256    expected contents of the three deposited archives
  FITS/, corpus/, corpus_item7/   populated by ./reproduce.sh fetch (or by refits)
Results/
  csv/, tex/, xlsx/    every table and in-text number as csv / LaTeX / workbook
  Figures/             figure trees (PDF + PNG), one per run label (Section 5)
  manuscript/          main-text figures + sidecar csv; stage3/ Stage 3 outputs
  manuscript_files.csv manuscript figure file -> package file, sha256, producer
  revision_rev2/       acceptance report of the revision batch (67 checks), inputs, manifest
```

## 2. Requirements

**Software.** R 4.6.1 (R >= 4.6 required). `Rscript install.R` installs the CRAN packages from the Posit Package Manager snapshot of 2026-07-01, which serves exactly the versions that produced the shipped results (all 70 packages of the dependency closure were checked; list in `sessionInfo.txt`), and OpTop 0.20.1 (commit `cba1273ea8`) and NLPstudio 1.2.0 (commit `771b42e328`) from GitHub; installed copies are accepted only at those versions and commits. `Rscript install.R lib=.Rlib` installs into a fresh project-local library, which `reproduce.sh` then uses automatically; this is the recommended route because packages already installed elsewhere are not downgraded.

OpTop compiles from source (C++17 with RcppArmadillo and OpenMP): it needs a C++ toolchain and a Fortran runtime (macOS: Xcode Command Line Tools and the CRAN gfortran build; Linux: build-essential and gfortran; Windows: Rtools). The revision batch runner needs Python >= 3.8 (standard library only). `reproduce.sh` needs bash, `curl` and `shasum`.

**Troubleshooting a source installation.** If a conda or miniconda installation is on your `PATH`, building the `xml2` package from source can pick up conda's `libxml2` (through `pkg-config` / `xml2-config`), after which R cannot load it and `xml2`, `tm`, `topicmodels` and `quanteda` fail to install. Remove conda from `PATH` for the installation (for example `conda deactivate`, then make sure no `.../miniconda3/bin` entry remains) and rerun `Rscript install.R lib=.Rlib`; packages already installed are kept.

**Hardware.** Results were produced on an Apple M3 Pro (12 cores, 36 GB RAM), macOS 26.5 and 27.0; the revision batch used 8 workers.

**Disk.** Repository about 85 MB; result objects 272 MB; fit cache 2.7 GB; MD&A corpus 178 MB.

**Time** (8 workers on the hardware above): unit gates 10 s; `fetch results` a few seconds after download; `paper` about 1.5 minutes; `stage3` about 45 minutes (completion residuals 5 minutes, completion reference curves 40 minutes); `revision` about 28 hours; `smoke` under a minute; full refits of the simulations about a day and of the MD&A application about 7 hours.

## 3. Data

**Sources.** The empirical corpus consists of Item 7 (Management's Discussion and Analysis) sections of Form 10-K filings on the SEC's EDGAR system (public domain), fiscal years 2015-2016: a `quanteda` corpus of 14,817 filings filed 2016-17 with document variables `cik`, `fyear`, `sic`. `Code/prep_mdna.R` applies the paper's filters (13,309 fiscal-2015/16 filings; 13,120 after removing 189 duplicate firm-period observations; 11,491 after the 200-token length floor; split 8,025 / 3,466) and writes `Data/MDNA/mdna_prep_2015_2016.qs2`, which is shipped in this repository. Simulated data are generated by seeded code (`seed_base = 1970000`; every realised seed is in `Results/csv/seeds.csv`).

**Deposit.** Objects too large for the repository are deposited in a Zenodo data record (CC0-1.0) as three archives, fetched and verified by `./reproduce.sh fetch {results|fits|corpus|all}`. The expected sha256 of every member is tracked in `Data/MANIFEST-*.sha256`:

| Archive | Size | Contents | Needed by |
|---|---|---|---|
| `result_objects.tar` | 259 MB, 23 files | the 16 result objects of the revision batch (`*_rev2`; sha256 also in `Results/revision_rev2/output_manifest.csv`), the two rev1 rescoring objects, the completion reference objects (`e2_completion_truth_rev2*`), the Stage 3 moment matrices, the 300-document MD&A prep (revision smoke gate) | `paper`, `stage3`, `revision` |
| `fit_cache.tar` | 2.6 GB, 3,076 files | every fitted model (WarpLDA slim caches), identical to `Results/csv/fit_manifest_rev1.csv` | `stage3`, `revision`, unit-free reruns |
| `mdna_corpus.tar` | 178 MB | the raw Item 7 corpus `Data/corpus_item7/item7_2016_2017.qs2` | `full-mdna` |

`./reproduce.sh fetch WHAT --local DIR` takes the archives from a local directory instead of downloading them.

## 4. Lanes

| Lane | Command | Reads | Produces | Exact? |
|---|---|---|---|---|
| test | `./reproduce.sh test` | nothing | unit gates U1-U25 | — |
| **paper** (default) | `./reproduce.sh paper` | tracked caches + result objects | every exhibit (Section 6) | yes: csv/tex byte-identical, PDFs up to dates |
| stage3 | `./reproduce.sh stage3 8` | + fit cache | MD&A completion-residual tests (Table S25); per-fit completion reference curves (Table S5 completion block) | yes (fixed seeds; verified 9 Oct 2026) |
| revision | `./reproduce.sh revision 8` | + fit cache | the `*_rev2` result objects, acceptance report | yes with the deposited fits; completed stages whose inputs are unchanged are verified and skipped |
| legacy | `./reproduce.sh legacy` | tracked caches | the pre-revision exhibit trees (2 Sep 2026) | yes |
| smoke | `./reproduce.sh smoke` | nothing | tiny end-to-end run | — |
| full-sims / full-mdna / robustness | `./reproduce.sh full-sims 11` etc. | nothing / corpus | new fits and caches | statistically equivalent only |

**The `paper` lane** runs, in this order: `postprocess_rev2.R`, `postprocess_revision.R`, `make_revision_tex.R`, `make_revision_figures.R`, `postprocess_e4.R` and `postprocess_e2_gap.R` on the revision objects, `make_figures.R` for the revision tree (E1, E2, E4, E6; `in_suffix=_rev2`) and the pre-revision tree (E3, E5), `make_tables.R`, `make_mdna_outputs.R in_suffix=_rev2`, and the scripts of `Code/manuscript/` that need no fit cache. The overrides `K_true=40 W=5000 K_grid=10:100:10 fit_method=WarpLDA alpha=0.1 beta=0.01` identify the published simulation design (caches are keyed by design, not by label). When typing these commands by hand in zsh, write the overrides as separate words (zsh does not split a variable holding several of them).

**Refits.** WarpLDA is multithreaded collapsed Gibbs sampling, so refitted models are statistically equivalent but not bit-identical: fit-dependent third decimals can move and the revision batch's fit-cache gate (sha256 of all 3,076 fits) rejects a rebuilt cache by design. Exact reruns of the evaluation-only analyses therefore use the deposited fit cache.

**Revision batch.** `Code/run_revision_batch.sh <workers> {A|B|final|all} _rev2 [--dry-run]` runs stages 0-13 (smoke gate, MD&A design runs, E2, E1, E6, diagnostics; E4 both prior arms and the 10 x 10 selection replicates; post-processing and the 67-check acceptance report). Each stage records the hashes of the shared modules, `configs.R` and its own driver; a stage whose identity and outputs are intact is skipped. The smoke stage's identity also covers a legacy driver of the authors' working tree that is not shipped, so it reruns (about 45 s). The batch rewrites `Results/revision_rev2/inputs.json` when it starts.

## 5. Figure trees and output naming

| Tree under `Results/Figures/` | Content | Used by the manuscript |
|---|---|---|
| `warptestFULL_NULL_rev2/` | simulation figures from the revision objects (E1, E2, E4, E6) | Figures S1, S2, S6-S10 |
| `warptestFULL_NULL/` | simulation figures from the July caches (all studies) | Figures S3-S5 (E3), S11 (E5); the rest superseded |
| `mdna_2015_2016_rev2/` | MD&A figures from the revision object | Figure S12 |
| `mdna_2015_2016_rev1/`, `warptestFULL_NULL_rev1/` | revision phase rev1 | Figure S13 |
| `mdna_2015_2016/`, `dgp2/` | pre-revision MD&A figures; second generating configuration | not cited |
| `Results/manuscript/` | main-text Figures 1-4 | Figures 1-4 |

Files carrying `_rev1` / `_rev2` come from the September 2026 revision; files without a suffix come from the pre-revision runs (July and 31 Aug - 2 Sep 2026). Two scoping defects of the pre-revision runs were corrected in the revision, so where both exist **the suffixed file is the one the manuscript uses**: Test 3 strata at K = 40/50/60 (MD&A) and K\* = 40 (Study III) were computed with the strata of the smallest tested K, and held-out word-level indices used a null deviance without the Poisson linear term (Supplement Section S3). The pre-revision exhibits are kept for provenance and for the regression gates.

## 6. What produces what

Paths are relative to `Results/csv/` unless they start with `Results/`; `…` abbreviates the design tag `E*_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01` (Kstar20 tags: `E*_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01`). "Batch" means the result object written by the revision batch (Section 4). Producers are in `Code/`.

### 6.1 Main text

| Exhibit | Label | File / numerical source | Producer |
|---|---|---|---|
| Table 1 | `tab:tests_summary` | definitions of the three tests (no numbers) | `R/utils_moment_tests.R` |
| Figure 1 | `fig:e1_curves` | `Results/manuscript/F1_selection_support.pdf`; data `F1_seed_scores.csv`, `F1_panel_A.csv`, `F1_panel_B.csv` | `manuscript/make_main_figures.R` (E1 batch object; `rev2_unbinned_comparator_rev2.csv`) |
| Table 2 | `tab:e1_khat` | `rev2_selreps_selection_rules_rev2.csv`; `rev2_e1_matched_optimum_e1_base_rev2.csv`; perplexity `e1_selreps_khat_rev2_E1…csv`; NPMI `e1_khat_rev2_E1…csv`; Macro maximum `Results/manuscript/F1_seed_scores.csv`; Monte Carlo SEs `Results/manuscript/stage3/stage3_e1_table2_mcse.csv` | `postprocess_rev2.R`; `manuscript/make_main_figures.R`; `manuscript/stage3_e1_selreps.R` |
| Figure 2 | `fig:gap_channels` | `Results/manuscript/gap_decomposition.pdf`; data `gap_decomposition_data.csv`, `gap_decomposition_seed_data.csv` | `manuscript/make_main_figures.R` (from `e3_gap_decomposition_E3…csv`, `run_E3.R`) |
| Table 3 | `tab:e4_tests` | `rev2_e4_size_e4_base_rev2.csv`; `e4_power_sizeadj_rev2_E4…csv`; `rev2_e4_r2_e4_base_rev2.csv` | `postprocess_rev2.R`; `postprocess_e4.R` |
| Figure 3 | `fig:mdna_fit` | `Results/manuscript/F2_mdna_fit.pdf`; data `F2_data.csv`, `F2_gain_pairs.csv` (all 190 pairs), `F2_panel_B.csv` | `manuscript/make_main_figures.R` (MD&A batch object; `mdna_gain_profile_rev1.csv`) |
| Table 4 | `tab:mdna_results` | A: `mdna_selection_rules_rev1.csv`, `mdna_matched_optimum_rev1.csv`; B: `rev2_mdna_moment_tests_rev2.csv`, `rev2_mdna_residual_mass_rev2.csv` | `postprocess_revision.R`; `postprocess_rev2.R` |
| Figure 4 | `fig:mdna_mass` | `Results/manuscript/residual_mass.pdf`; data `residual_mass_data.csv`, `residual_contrast_data.csv`, `residual_cluster_sums.csv` (anonymous firm sums), `residual_interval_*.csv` | `manuscript/make_main_figures.R` |

### 6.2 Supplementary Material

| Exhibit | Label | File / numerical source | Producer |
|---|---|---|---|
| Table S1 | `tab:S_selection_full` | as Table 2, plus `rev2_e1_selection_rules_e1_base_rev2.csv` | `postprocess_rev2.R`; `manuscript/make_main_figures.R` |
| Table S2 | `tab:S_e1_r2` | `T1b_heldout_by_K_warptestFULL_NULL.csv` | `make_tables.R` |
| Figure S1 | `fig:S_e1_gains` | `Results/Figures/warptestFULL_NULL_rev2/E1/F2_adjacent_gains.pdf` | `make_figures.R` (revision tree) |
| Tables S3, S4 | `tab:S_e2_cal`, `tab:S_e2_cov` | `rev2_e2_coverage_e2_base_rev2.csv`, `rev2_e2_coverage_cells_e2_base_rev2.csv` | `postprocess_rev2.R` |
| Figure S2 | `fig:S_e2_cov` | `Results/Figures/warptestFULL_NULL_rev2/E2/F3_coverage_qq.pdf` | `make_figures.R` (revision tree) |
| Table S5 | `tab:S_e2_falsecert` | `Results/manuscript/stage3/`: `stage3_falsecert_fits.csv`, `stage3_falsecert_table.csv`, `stage3_e2_reference_curves.csv`, `stage3_completion_reference_curves.csv`, `stage3_completion_falsecert.csv`, `stage3_completion_paired_se.csv` (notes); LaTeX rows `tabS5_rows.tex` (reconstruction panels) | `manuscript/stage3_e2_falsecert.R`, `stage3_e1_selreps.R`, `stage3_make_tables.R`; completion block `stage3_completion_reference.R` (fit cache), `stage3_completion_paired_se.R` |
| Table S6 | `tab:S_unbinned` | `rev2_unbinned_comparator_rev2.csv`, `…_optima_rev2.csv`, `…_floor_sensitivity_rev2.csv` | `postprocess_rev2.R` |
| Table S7 | `tab:S_gap_channels` | `e3_gap_decomposition_E3…csv` (ten-seed means in `Results/manuscript/gap_decomposition_data.csv`) | `run_E3.R`; `manuscript/make_main_figures.R` |
| Tables S8, S9 | `tab:S_e3_gap_main`, `tab:S_e3_len` | `T6a_gap_decomposition_…`, `T6b_length_stats_warptestFULL_NULL.csv` | `make_tables.R` |
| Figures S3-S5 | `fig:S_e3_gap_curves`, `fig:S_e3_gap_decomp_exact`, `fig:S_e3_scatter` | `Results/Figures/warptestFULL_NULL/E3/F6_gap_ci.pdf`, `F6b_gap_decomposition.pdf`, `F7_doc_scatters.pdf` | `make_figures.R` (E3) |
| Table S10 | `tab:S_e6` | `rev2_e6_word_curves_rev2.csv`, `rev2_e6_lemma_crosspath_rev2.csv` | `postprocess_rev2.R` |
| Figures S6, S7 | `fig:S_e6_gap`, `fig:S_e6_freq` | `Results/Figures/warptestFULL_NULL_rev2/E6/F11_word_micro_macro.pdf`, `F12_word_fit_vs_freq.pdf` | `make_figures.R` (revision tree) |
| Table S11 | `tab:S_e4_size` | `rev2_e4_size_e4_base_rev2.csv`, `rev2_e4_size_e4_prior_rev2.csv` (`e4_size_rev2_E4…_fa0p1/fa0p5…csv`) | batch (E4, both prior arms); `postprocess_rev2.R` |
| Table S12 | `tab:S_e4_power` | `rev2_e4_power_e4_base_rev2.csv`, `rev2_e4_r2_e4_base_rev2.csv` | `postprocess_rev2.R` |
| Table S13 | `tab:S_e4_sizeadj` | `e4_power_sizeadj_rev2_E4…csv` (`power_adj_seed`) | `postprocess_e4.R` |
| Figures S8, S9 | `fig:S_e4_power`, `fig:S_e4_fitvstests` | `Results/Figures/warptestFULL_NULL_rev2/E4/F4_power_curves.pdf`, `F5_fit_vs_tests.pdf` | `make_figures.R` (revision tree) |
| Table S14 | `tab:S_e4_words` | row 1 `T5_planted_words_warptestFULL_NULL.csv`; row 2 `rev2_e4_planted_words_e4_base_rev2.csv` | `make_tables.R`; `postprocess_rev2.R` |
| Figure S10 | `fig:S_e4_words_fig` | `Results/Figures/warptestFULL_NULL_rev2/E4/F5b_word_ranks.pdf` | `make_figures.R` (revision tree) |
| Table S15 | `tab:S_e5` | `T7a_grid_sensitivity_…`, `T7b_minbin_warptestFULL_NULL.csv` | `make_tables.R` |
| Figure S11 | `fig:S_e5_sens` | `Results/Figures/warptestFULL_NULL/E5/F9_design_sensitivity.pdf` | `make_figures.R` (E5) |
| Table S16 | `tab:S_resolution` | `rev2_support_resolution_{e1_base,e1_dgp2,e5,mdna}_rev2.csv` | `postprocess_rev2.R` |
| Table S17 | `tab:S_dgp2` | `rev2_e1_selection_rules_e1_dgp2_rev2.csv`; NPMI/perplexity rows `e1_khat_rev2_E1_full_Kstar20…csv` | `postprocess_rev2.R` |
| Tables S18, S19 | `tab:S_mdna_gains`, `tab:S_mdna_sel_family` | `mdna_gain_profile_rev1.csv`, `mdna_selection_rules_rev1.csv` (`Results/tex/rev1_S13_gain_profile.tex`, `rev1_S13b_selection_by_family.tex`) | `postprocess_revision.R`; `make_revision_tex.R` |
| Figure S12 | `fig:S_mdna_gap` | `Results/Figures/mdna_2015_2016_rev2/MDNA/R1b_micro_macro_gap.pdf` | `make_mdna_outputs.R in_suffix=_rev2` |
| Table S20 | `tab:S_mdna_word` | `mdna_wordlevel_conventions_rev1.csv`; gaps `R2_battery_mdna_2015_2016.csv`, `Results/manuscript/stage3/stage3_mdna_gap_cluster.csv` | `rescore_revision.R` (fit cache); `make_mdna_outputs.R`; `manuscript/stage3_mdna_subsamples.R` |
| Table S21 | `tab:S_mdna_wcurve` | `rev2_mdna_word_curves_rev2.csv` | `postprocess_rev2.R` |
| Figure S13 | `fig:S_mdna_tests` | `Results/Figures/mdna_2015_2016_rev1/MDNA/R4_moment_diagnostics.pdf` | `make_revision_figures.R` |
| Tables S22-S24 | `tab:S_mdna_joint`, `tab:S_mdna_strata`, `tab:S_mdna_mass` | `mdna_moment_tests_rev1.csv`, `mdna_moment_tests_cluster_rev1.csv`, `mdna_moment_strata_rev1.csv`, `mdna_residual_mass_rev1.csv` (revision twins `rev2_mdna_*_rev2.csv`); LaTeX bodies `Results/tex/rev1_S16*.tex` | `rescore_revision.R`; `make_revision_tex.R` |
| Table S25 | `tab:S_mdna_completion` | `Results/manuscript/stage3/stage3_mdna_completion_tests.csv` (+ `_strata`, `_masses`); rows `tabS25_rows.tex` | `manuscript/stage3_mdna_completion_moments.R` (fit cache); `stage3_make_tables.R` |
| Table S26 | `tab:S_mdna_vocab` | `Results/tex/rev1_S17_vocabulary.tex` (`mdna_wordlevel_worst_fit_paper_null_rev1.csv`) | `make_revision_tex.R` |
| Table S27 | `tab:S_mdna_boiler` | `A2_boilerplate_by_industry_mdna_2015_2016.csv` | `make_mdna_outputs.R` |
| Table S28 | `tab:S_mdna_sens` | `rev2_mdna_design_sensitivity_rev2.csv`; `mdna_selection_firm_sensitivity_rev1.csv`; present-firm rows `Results/manuscript/stage3/stage3_mdna_subsample_selections.csv` (`present_rows.tex`) | `postprocess_rev2.R`; `postprocess_revision.R`; `manuscript/stage3_mdna_subsamples.R` |
| Table S29 | `tab:S_mdna_delta` | `mdna_delta_sensitivity_rev1.csv` | `postprocess_revision.R` |

### 6.3 Numbers in the text

Every number in the text is a row of one of the files above. Additional sources:
- Section 5.1 and S5.1 (per-fit adequacy, false certifications, fit 7 and fit 8): `Results/manuscript/stage3/stage3_e2_reference_curves.csv`, `stage3_completion_reference_curves.csv`, `stage3_completion_paired_se.csv`.
- Sections 6.2-6.3 and S6.2-S6.6 (firm-clustered gap intervals, present/absent firms, bootstrap interval for the displaced mass H, completion residuals, the decomposition W = W_shift + W_homogeneity): `Results/manuscript/stage3/stage3_mdna_gap_cluster.csv`, `stage3_mdna_subsample_*.csv`, `stage3_mdna_H_bootstrap.csv`, `stage3_mdna_completion_*.csv`, `stage3_mdna_s64_wald_subtests.csv`.
- Section 6.1 (restart dispersion): `rev2_mdna_restart_range_rev2.csv`, `rev2_mdna_restart_paired_rev2.csv`.
- Coverage and selection by evaluation size: `rev2_e2_coverage_e2_base_rev2.csv`, `rev2_e2_selection_by_Jev_e2_base_rev2.csv`; second configuration `rev2_e2_*_e2_dgp2_rev2.csv`.

`Results/manuscript/manuscript_checks.json` records the 524 machine checks of the final LaTeX sources against these files (9 October 2026; all passed), with the sha256 of every file checked.

## 7. Versions and provenance

The cached results were produced in two environments (details and per-object stamps in `sessionInfo.txt`): the July 2026 baseline caches with OpTop 0.14.0/0.14.1 and NLPstudio 1.1.1, and everything from 31 August 2026 onwards in one CRAN environment with OpTop 0.20.0 (robustness arms, rev1 rescoring) or 0.20.1 (revision batch, Stage 3, final regeneration). OpTop's changelog records two bit-level numerical changes between 0.14.1 and 0.20.0 (at most about 1e-12 relative); `Code/R/utils_heldout.R` rebuilds the dense partition from OpTop >= 0.15's compressed form, and the unit gates enforce agreement with OpTop's compiled path (U1/U2 at 1e-10, U15 at 1e-8). OpTop 0.20.1 adds only a classed warning when its word-level null is used with an external baseline; the replication code computes that null itself (Supplement Section S3) and gate U20 checks the warning, so the gates require 0.20.1.

## 8. Verification record

On 9 October 2026, with the environment of `install.R` (R 4.6.1, macOS 27.0.1), a clean copy of this repository was checked using only the commands above:
- `./reproduce.sh paper` refused to run before the result objects were present; `./reproduce.sh fetch results|fits|corpus --local` then unpacked the three archives and verified every member against `Data/MANIFEST-*.sha256`.
- `./reproduce.sh test`: all unit gates (U1-U25) pass.
- `./reproduce.sh paper` (83 s) followed by the Stage 3 scripts (`stage3_mdna_completion_moments.R`, whose internal gate reproduces the stored reconstruction statistics to 1e-6; the completion reference curves in smoke mode; `stage3_completion_paired_se.R`; `stage3_make_tables.R`): of the 330 files under `Results/`, 300 were regenerated byte for byte and the other 30 differ only in metadata, namely 28 PDFs (24 identical after removing creation dates and document IDs, 4 identical when rasterised) and the 2 Excel workbooks (identical sheets). The 40-minute full completion-reference run was checked through its smoke and pilot modes (byte-identical) and by recomputing every curve from the saved per-document scores (agreement to 6e-16).
- `bash Code/run_revision_batch.sh 8 all _rev2 --dry-run` with the deposited fit cache reports every stage of the revision batch as valid except the smoke gate (Section 4).
- `./reproduce.sh smoke` completes.
- `Rscript install.R lib=.Rlib` in a fresh library (only R's base packages visible) installed all 70 dependencies at the snapshot versions and both GitHub packages at their commits; with that library `./reproduce.sh test` passed and `./reproduce.sh paper` regenerated the same 300 files byte for byte, PNG previews included (they are drawn by `ragg`, pinned in `install.R`). On this machine the source build of `xml2` needed conda removed from `PATH` (Section 2).
- Every figure file of the manuscript is in this package, byte-identical to the submitted file (16 of 17) or identical up to PDF metadata (Figure S10): `Results/manuscript_files.csv`.
- 524 machine checks of the manuscript's numbers against this package pass (`Results/manuscript/manuscript_checks.json`).

## 9. Notes for referees

- `Code/README.md` is the development manual of the pipeline (CLI grammar, override keys, engines, monitoring, design notes); where its examples use exploratory profiles, the commands of **this** README and of `reproduce.sh` are the ones behind the manuscript.
- Unit gates (`Code/tests_unit.R`): U1-U2 partition and scoring equal OpTop's compiled path; U3-U4 gap identity and channels (Proposition 1); U5 token split; U6 selection rule; U7 instruments and Wald size; U8 fold-in; U9 WarpLDA cache round trip; U10 Lemma S1 across code paths; U11 word-level helpers; U12-U13 prior handling; U14 fit-cache keys; U15 native partition and held-out scorer; U16 null-discrepancy floor; U17 K-specific Test 3 instruments; U18 total-gain and simultaneous selectors, "none certified", cluster SEs; U19 Poisson-form held-out word null; U20 OpTop 0.20.1 warning contained; U21 threshold c and floor delta independent; U22 checkpoint identity; U23 word-null convention; U24 support resolution; U25 completion scores only scored tokens.
- The null-discrepancy floor (`min_null`, Remark 1) is applied uniformly; excluded shares are reported (`null_excl_share`, MD&A boilerplate tables).

## 10. Citation and licences

Please cite the paper and this package (`CITATION.cff`; GitHub's "Cite this repository"). Code: MIT. Data, result objects, tables and figures: CC0-1.0 (`LICENSE-DATA`). SEC filing text is in the public domain.
