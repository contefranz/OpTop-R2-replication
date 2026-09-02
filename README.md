# Replication Package

**Goodness-of-Fit Indices and Diagnostics for Topic Models**
Craig M. Lewis (Vanderbilt University) and Francesco Grossetti (Bocconi University)
Submitted to the *Journal of the Royal Statistical Society, Series B*. Contact: francesco.grossetti@unibocconi.it
Repository: https://github.com/contefranz/OpTop-R2-replication

This package reproduces **every table and figure** in the main paper and in the online Supplementary Material: the four simulation studies (Section 5) and the MD&A empirical application (Section 6). All methods are implemented in the R package [`OpTop`](https://github.com/contefranz/OpTop); this package contains the experiment drivers, the empirical pipeline, the cached results behind the published exhibits, and the data.

Exhibit numbering below refers to the September 2026 manuscript sources (the JRSS-B submission); every entry also carries its LaTeX label, so the mapping survives renumbering.

---

## 1. Contents

```
.here                 project-root anchor (do not delete: paths resolve against it)
README.md             this file
LICENSE               MIT licence for the code (SEC filing text is public domain)
CITATION.cff          citation metadata (GitHub "Cite this repository")
.zenodo.json          metadata for the Zenodo archive of each GitHub release
install.R             installs the pinned package versions (run once)
reproduce.sh          one-command entry points for the lanes below
sessionInfo.txt       the two environments behind the shipped results (July 2026
                      baseline, September 2026 arms) + per-cache version stamps
Code/
  R/                  analysis modules (I/O, DGP, fitting, held-out scoring,
                      inference, moment tests, comparators, plot theme)
  config/configs.R    all design parameters, profiles (smoke/pilot/full), seeds
  run_E1.R..run_E6.R  simulation drivers (idempotent, cache-backed)
  make_all.R          orchestrator: E1-E6 + figures + tables
  make_figures.R      figures from saved results only
  make_tables.R       tables (CSV + LaTeX + one Excel workbook)
  prep_mdna.R         MD&A corpus -> analysis-ready document-term objects
  run_mdna.R          MD&A estimation + all analyses (fit, battery, tests);
                      c_part=/out_suffix= support the c- and grid-sensitivity
                      re-evaluations without touching the baseline results
  make_mdna_outputs.R MD&A paper figures and tables
  run_E1_selreps.R    Study-I selection on S x R_sel replicated evaluation
                      corpora (evaluation-only; cached fits)
  run_mdna_restarts.R across-restart dispersion of the MD&A held-out index
  postprocess_e4.R    size-adjusted power (empirical critical values) +
                      contamination-refit seed profile, from saved E4 output
  postprocess_e2_gap.R gap-CI coverage diagnosis, from saved E2 output
  run_tonight.sh      one-shot robustness batch (all of the above + second-DGP
                      and correct-prior arms), pilot-gated and logged; stages
                      2-7 are resilient (a failure is logged, the chain
                      continues) and `run_tonight.sh <workers> <from>` resumes
                      from a given stage id
  list_runs.R         lists every cached run + the exact command to re-plot it
  tests_unit.R        correctness gates (non-zero exit on failure)
  README.md           full pipeline manual (CLI grammar, override keys,
                      engines, monitoring, design notes)
Data/
  E1/..E6/            simulation result caches: the definitive July 2026 run
                      [14 MB] plus the September 2026 robustness arms:
                      replicated selection (E1), second DGP (E1, E2),
                      matched-prior size arm (E4) [4 MB]
  MDNA/               mdna_prep_2015_2016.qs2   analysis-ready corpus [20 MB] (shipped; also DEPOSIT)
                      mdna_results_MDNA_2015_2016.qs2  all MD&A results [9 MB]
                      mdna_results_MDNA_2015_2016_{c05,c2,g100}.qs2
                                                c/grid sensitivity re-evaluations [28 MB]
  corpus_item7/       item7_2016_2017.qs2  raw Item-7 corpus [178 MB] (DEPOSIT)
  FITS/, corpus/      empty; populated only by full recomputes
Results/
  csv/, tex/, xlsx/   the published tables (+ seeds.csv: every seed used)
  Figures/            the published figures (PDF + PNG); dgp2/ holds the
                      second-DGP analogues of Figures 1, S1, S2 (not cited)
  FIGURE_GUIDE.md     plain-language guide to reading each simulation figure
```

Files marked **(DEPOSIT)** are data, not code: they are deposited on Zenodo with their own DOI (see Section 7). The raw corpus exceeds GitHub's 100 MB file limit and is excluded from the GitHub repository by `.gitignore`; the analysis-ready object is small enough to ship here as well, so Lanes B and C run from a clone alone.

## 2. Requirements

R ≥ 4.6 with the packages pinned in `install.R` / `sessionInfo.txt`. The two non-CRAN dependencies are pinned to the commits of the environment that produced the September 2026 robustness arms and under which every shipped exhibit was last regenerated (2 September 2026):

```r
remotes::install_github("contefranz/OpTop@1166daed33")      # v0.20.0 (main, 2026-07-22)
remotes::install_github("contefranz/NLPstudio@771b42e328")  # v1.2.0  (main, 2026-07-23)
```

**Provenance.** The shipped caches were produced in two environments; `sessionInfo.txt` records both, together with the version stamp each cache carries in `attr(x, "run_meta")$packages`:

* the baseline simulation caches E1–E6 and the MD&A baseline — 11–14 July 2026, OpTop 0.14.1 (E4 and E6: 0.14.0) with NLPstudio 1.1.1;
* the September 2026 robustness arms (replicated selection, second DGP, matched-prior size arm, MD&A c/grid sensitivity, restart dispersion), all post-processing, and every shipped figure and table — OpTop 0.20.0 with NLPstudio 1.2.0.

**Version compatibility.** The pipeline runs unchanged under both environments. OpTop ≥ 0.15 returns the harmonised partition in compressed form; `Code/R/utils_heldout.R` rebuilds the dense mask transparently (`.densify_partition()`), and the unit gates (`Code/tests_unit.R`, U1/U2/U15/U16) enforce agreement with the frozen pure-R reference implementation and with the package's compiled path at 1e-10. Between 0.14.1 and 0.20.0 OpTop's changelog records two bit-level numerical changes (at most about 1e-12 relative; releases 0.15.0 and 0.17.0), so every shipped number is reproduced to that tolerance under either version: the gates pass under 0.20.0, and Lane A below rebuilds the shipped exhibits from the caches under 0.20.0. The null-discrepancy floor introduced in 0.14.1 (`min_null`, Remark 1) does not touch the two 0.14.0 caches: E4 and E6 score under held-out reconstruction and in-sample, where no document falls below the floor (it binds only under completion), and the word-level indices are unaffected by construction. NLPstudio 1.2.0 is a performance release with identical results; it enters the pipeline only through the `fit_topic_model()` facade.

Shipped results were produced with R 4.6.1 on an Apple M3 Pro (12 cores, 36 GB RAM), macOS 26.5.2 (July) and 27.0 (September). Run `Rscript install.R` once; it checks versions and installs anything missing.

**Always run commands from this directory** (the folder containing `.here`).

## 3. Verify first

```sh
Rscript Code/tests_unit.R        # unit gates: partition/index/fold-in/moment-test
                                 # correctness, incl. equality against OpTop's
                                 # compiled implementations; non-zero exit on failure
```

## 4. Lane A — regenerate every exhibit from the shipped caches (~minutes)

No model is refit; figures and tables are rebuilt from `Data/*/…qs2` and overwrite `Results/` in place.

```sh
./reproduce.sh exhibits
```

which runs (with the script's `SIM_DESIGN` variable expanded):

```sh
# Simulations (Section 5 + Supplement S5): figures then tables
Rscript Code/make_figures.R full K_true=40 W=5000 K_grid=10:100:10 \
        fit_method=WarpLDA alpha=0.1 beta=0.01 label=warptestFULL_NULL
Rscript Code/make_tables.R  full K_true=40 W=5000 K_grid=10:100:10 \
        fit_method=WarpLDA alpha=0.1 beta=0.01 label=warptestFULL_NULL

# Post-processing of the saved E4 and E2 objects (Table 2 parentheses,
# Table S9, the gap-CI diagnosis of S5.1 and its second-DGP check)
Rscript Code/postprocess_e4.R     file=Data/E4/e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2
Rscript Code/postprocess_e2_gap.R file=Data/E2/e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2
Rscript Code/postprocess_e2_gap.R file=Data/E2/e2_results_E2_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01.qs2

# Second-DGP figures (Supplement S5.5; produced, not cited)
Rscript Code/make_figures.R full K_true=20 J_train=1000 W=10000 L=500 K_grid=5:50:5 \
        fit_method=WarpLDA alpha=0.1 beta=0.01 exp=E1,E2 label=dgp2

# Empirical application (Section 6)
Rscript Code/make_mdna_outputs.R y1=2015 y2=2016
```

The `key=value` overrides identify *which cached run* to load (caches are keyed by design parameters, never by label); `label=` sets *where* outputs go. `Rscript Code/list_runs.R` prints the same commands by scanning the caches.

## 5. Lane B — full recompute from scratch

### 5.1 Simulations (budget ~1 day of wall-clock on 11 workers)

```sh
Rscript Code/make_all.R full 11 \
    K_true=40 W=5000 K_grid=10:100:10 fit_method=WarpLDA alpha=0.1 beta=0.01 \
    J_eval=5000 S=10 R_eval=10 J_truth=2000 J_ev=250 R_power=250 J_center=10000 K_test=10,40 \
    label=warptestFULL_NULL
```

This is the definitive published design, reconstructed from the configuration embedded in the shipped caches: generating LDA with $K^\star=40$, $W=5{,}000$, $J_{tr}=1{,}000$, $\alpha_{DGP}=0.5$, $\beta_{DGP}=0.01$; WarpLDA estimation (3 restarts) with deliberately misspecified fitted priors ($\alpha=0.1$, $\beta=0.01$); grid $K=10{:}100$ by 10; ten seeds. Shared overrides apply where each design uses them and are ignored elsewhere (`S` sizes E1/E3/E6; `S_train`/`R_eval`/`R_null`/`R_power` size the conditional Monte Carlo of E2/E4; `J_ev`/`J_center`/`K_test` are E4's). Seeds derive deterministically from `seed_base=1970000` (see `Code/config/configs.R`; every realized seed is logged to `Results/csv/seeds.csv`).

Cost profile: LDA fitting is the small share (Table `T8`: 33 s/fit average under WarpLDA; the E1 fit pool is ~55 core-minutes); the held-out scoring passes dominate. Progress streams to `Data/FITS/prefit.log` and `Data/scoring.log`. Every stage is idempotent — interrupted runs resume from the caches.

A fast end-to-end exercise of the identical code paths (tiny corpus, ~2 minutes):

```sh
./reproduce.sh smoke        # = Rscript Code/make_all.R smoke
```

### 5.2 Empirical application (~7 hours on 6 workers)

```sh
Rscript Code/prep_mdna.R y1=2015 y2=2016 workers=6                       # ~minutes
Rscript Code/run_mdna.R  y1=2015 y2=2016 workers=6 K_grid=10:200:10 refine_span=0   # ~6.5 h
Rscript Code/make_mdna_outputs.R y1=2015 y2=2016                          # <1 min
```

`prep_mdna.R` starts from the raw corpus `Data/corpus_item7/item7_2016_2017.qs2` (fiscal-year filter, one filing per firm-period, tokenization, training-vocabulary pruning at document frequency ≥ 5, length floor 200, seeded 70/30 train/evaluation split) and writes `Data/MDNA/mdna_prep_2015_2016.qs2` — the shipped copy of which lets you skip this step. The driver fits WarpLDA (3 restarts, engine-default priors) on $K=10{:}200$ by 10 and runs the full battery; moment tests run at $\hat K$ and its grid neighbors (`tests_K=hat`; pass `tests_K=all` to test every grid point, roughly one extra scoring pass).

**Determinism.** All randomness is seeded, but WarpLDA is multithreaded collapsed Gibbs: refit results are statistically equivalent, not bit-identical (fit-dependent third decimals can move; selections and test conclusions are stable across restarts). Lane A is exactly reproducible; the shipped caches are the published numbers.

### 5.3 Robustness arms (Lane C; about 16 h on 8 workers)

```sh
./reproduce.sh robustness 8        # = bash Code/run_tonight.sh 8
```

`Code/run_tonight.sh` chains, cheapest first, the September 2026 arms behind the $n=100$ rows of Table 1, the lower panel of Table S7, Tables S12 and S19, and the in-text sentences of Sections 5.2, 5.3, 6.1 and 6.2: the unit gates and a pilot gate for the second DGP (stages 0–1, fatal on failure); the MD&A restart dispersion at $K=50$ (`run_mdna_restarts.R`); the matched-prior E4 size arm (`run_E4.R … alpha=0.5 S_train=10 R_null=500 R_power=0 alternatives=none`); the replicated Study-I selection (`run_E1_selreps.R … S=10 R_sel=10 J_eval=5000`, evaluation-only on the cached fits); the second generating configuration (`run_E1.R` / `run_E2.R` with `K_true=20 J_train=1000 W=10000 L=500 K_grid=5:50:5`); the MD&A c- and grid-sensitivity re-evaluations (`run_mdna.R … c_part=0.5 out_suffix=_c05`, `c_part=2 out_suffix=_c2`, `K_grid=10:100:10 out_suffix=_g100`); and finally the post-processing and figure/table regeneration of Lane A (stage 7). Stages 2–7 are resilient (a failure is logged and the chain continues) and `bash Code/run_tonight.sh <workers> <from-stage>` resumes from a stage id. The outputs are the September caches and csv files listed in Section 6; the shipped copies were produced by this lane on 1–2 September 2026 (16 h wall-clock on 8 workers) under OpTop 0.20.0 / NLPstudio 1.2.0.

## 6. What produces what

### 6.1 Main text

| Exhibit | LaTeX label | Output file | Produced by | Computed by |
|---|---|---|---|---|
| Figure 1 | `fig:e1_curves` | `Results/Figures/warptestFULL_NULL/E1/F1_heldout_curves.pdf` | `make_figures.R` | `run_E1.R` |
| Figure 2 | `fig:mdna_fit` | `Results/Figures/mdna_2015_2016/MDNA/R1_fit_curves.pdf` | `make_mdna_outputs.R` | `run_mdna.R` |
| Figure 3 | `fig:mdna_tests` | `Results/Figures/mdna_2015_2016/MDNA/R4_moment_diagnostics.pdf` | `make_mdna_outputs.R` | `run_mdna.R` |
| Table 1 | `tab:e1_khat` | n=100 rows: `Results/csv/e1_selreps_khat_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.csv`; n=10 rows: `Results/csv/T1_khat_distribution_warptestFULL_NULL.csv` | `run_E1_selreps.R` / `make_tables.R` | `run_E1_selreps.R` / `run_E1.R` |
| Table 2 | `tab:e4_tests` | `T3_moment_size_…csv` + `T4_moment_power_…csv` + `T4b_r2_under_misspec_…csv`; size-adjusted values in parentheses: `e4_power_sizeadj_…csv` (`power_adj_seed` at K=40) | `make_tables.R` + `postprocess_e4.R` | `run_E4.R` |
| Table 3 | `tab:mdna_sel` | `Results/csv/R1_selection_mdna_2015_2016.csv` | `make_mdna_outputs.R` | `run_mdna.R` |
| Table 4 | `tab:mdna_tests` | `Results/csv/R2_battery_mdna_2015_2016.csv` (rows `block == "moment_tests"`); the full over-$K$ table including the $K=40/60$ columns is the `tests_ho` element of `Data/MDNA/mdna_results_MDNA_2015_2016.qs2` (drawn in the moment-diagnostics figure, right panel) | `make_mdna_outputs.R` | `run_mdna.R` |

The September 2026 draft moved the MD&A gap figure and the residual-vocabulary table to Supplement S6 and completed S6 with new exhibits. All are drawn from `Data/MDNA/mdna_results_MDNA_2015_2016.qs2` (elements named below) or the listed csv, producers `run_mdna.R` (+ `make_mdna_outputs.R` for the figure): `fig:S_mdna_gap` = `MDNA/R1b_micro_macro_gap.pdf`; Table S13 `tab:S_mdna_gains` = `gains`; Table S14 `tab:S_mdna_word` = `wstar`; Table S15 `tab:S_mdna_wcurve` = `wcurve`; Table S16 `tab:S_mdna_strata` = `strata_ho`; Table S17 `tab:S_mdna_vocab` = `Results/csv/R3_vocabulary_mdna_2015_2016.csv`; Table S18 `tab:S_mdna_boiler` = `Results/csv/A2_boilerplate_by_industry_mdna_2015_2016.csv`; Table S19 `tab:S_mdna_sens` = the `c05`/`c2`/`g100` results objects, selections recomputed on the coarse grid.

In-text numbers of Section 6 (fit levels, $\hat K$ by rule, CI bounds, gap channels, boilerplate shares, word-dual summaries, test statistics) are all rows of `R1_selection`/`R2_battery`/`R3_vocabulary`; the same content is browsable in `Results/xlsx/MDNA_tables_2015_2016.xlsx`, one sheet per table. Three September 2026 in-text additions: the Section 6.1 restart-dispersion sentence (< 0.005) = `Results/csv/mdna_restart_dispersion_MDNA_2015_2016.csv` (`run_mdna_restarts.R`); the Section 6.2 c/grid-sensitivity sentence = `Data/MDNA/mdna_results_MDNA_2015_2016_{c05,c2,g100}.qs2` (`run_mdna.R` with `c_part=`/`out_suffix=`/`K_grid=`, selections recomputed on the coarse grid); the Section 5.3 matched-prior sentence = `Results/csv/e4_size_…_fa0p5_fb0p01.csv` (`run_E4.R … alpha=0.5 R_null=500 R_power=0 alternatives=none`).

### 6.2 Supplementary Material (Section S5; Study I = E1–E2, Study II = E3+E6, Study III = E4, Study IV = E5)

| Exhibit | LaTeX label | Output file(s) (`Results/csv/…_warptestFULL_NULL.csv`, figures under `Results/Figures/warptestFULL_NULL/`) | Computed by |
|---|---|---|---|
| Table S1 | `tab:S_tests_summary` | static summary of the three tests — no computation | — |
| Table S2 | `tab:S_e1_r2` | `T1b_heldout_by_K` | `run_E1.R` |
| Table S3 | `tab:S_e2_cov` | `T2a_coverage` (notes draw on `T2b_gain_test_size`, `T2c_gain_test_power`, `T2e_khat_by_Jev`; `T2d_gap_coverage` backs the gap-CI claim) | `run_E2.R` |
| Table S4 | `tab:S_e3_gap_main` | `T6a_gap_decomposition` | `run_E3.R` |
| Table S5 | `tab:S_e3_len` | `T6b_length_stats` | `run_E3.R` |
| Table S6 | `tab:S_e6` | `T9_word_micro_macro` | `run_E6.R` |
| Table S7 | `tab:S_e4_size` | `T3_moment_size` (upper panel); lower matched-prior panel: `e4_size_…_fa0p5_fb0p01.csv` | `run_E4.R` (twice: baseline and `alpha=0.5` arm) |
| Table S8 | `tab:S_e4_power` | `T4_moment_power` + `T4b_r2_under_misspec` | `run_E4.R` |
| Table S9 | `tab:S_e4_sizeadj` | `e4_power_sizeadj_…_fa0p1_fb0p01.csv` (`power_adj_seed`, K=40) | `postprocess_e4.R` (from `run_E4.R` output) |
| Table S10 | `tab:S_e4_words` | `T5_planted_words` (the per-stratum illustration its notes cite = `T5b_strata_illustration`) | `run_E4.R` |
| Table S11 | `tab:S_e5` | `T7a_grid_sensitivity` + `T7b_minbin` | `run_E5.R` |
| Table S12 | `tab:S_dgp2` | `e1_khat_E1_full_Kstar20_J1000_W10000_…csv` (+ coverage sentence: `e2_coverage_…Kstar20…csv`, `e2_gap_diagnosis_…Kstar20…csv`, `e2_khat_distribution_…Kstar20…csv`) | `run_E1.R` / `run_E2.R` with the `dgp2` overrides of `run_tonight.sh` |
| Figure S1 | `fig:S_e1_gains` | `E1/F2_adjacent_gains.pdf` | `run_E1.R` |
| Figure S2 | `fig:S_e2_cov` | `E2/F3_coverage_qq.pdf` | `run_E2.R` |
| Figure S3 | `fig:S_e3_gap_curves` | `E3/F6_gap_ci.pdf` | `run_E3.R` |
| Figure S4 | `fig:S_e3_gap_decomp_exact` | `E3/F6b_gap_decomposition.pdf` | `run_E3.R` |
| Figure S5 | `fig:S_e3_scatter` | `E3/F7_doc_scatters.pdf` | `run_E3.R` |
| Figure S6 | `fig:S_e6_gap` | `E6/F11_word_micro_macro.pdf` | `run_E6.R` |
| Figure S7 | `fig:S_e6_freq` | `E6/F12_word_fit_vs_freq.pdf` | `run_E6.R` |
| Figure S8 | `fig:S_e4_power` | `E4/F4_power_curves.pdf` | `run_E4.R` |
| Figure S9 | `fig:S_e4_fitvstests` | `E4/F5_fit_vs_tests.pdf` | `run_E4.R` |
| Figure S10 | `fig:S_e4_words_fig` | `E4/F5b_word_ranks.pdf` | `run_E4.R` |
| Figure S11 | `fig:S_e5_sens` | `E5/F9_design_sensitivity.pdf` | `run_E5.R` |

All figures are also produced as PNG next to the PDFs. Simulation tables ship additionally as LaTeX fragments (`Results/tex/`) and as one workbook (`Results/xlsx/Section5_tables_warptestFULL_NULL.xlsx`).

### 6.3 Produced and shipped, not (yet) cited

`MDNA/R1c_adjacent_gains.pdf` (figure version of Table S13), `MDNA/A1_word_micro_macro.pdf` (figure version of Table S15), `MDNA/A3_fit_vs_length.pdf` (document fit vs length), and `E1/F1b_micro_curves.pdf` (Micro analogue of Figure 1) are produced but not cited. So is the second-DGP figure tree `Results/Figures/dgp2/` (`E1/F1_heldout_curves.pdf`, `E1/F1b_micro_curves.pdf`, `E1/F2_adjacent_gains.pdf`, `E2/F3_coverage_qq.pdf`: the analogues of Figures 1, S1 and S2 for the configuration of Table S12).

## 7. Data

**Sources.** The empirical corpus consists of Item 7 (Management's Discussion and Analysis) sections extracted from Form 10-K filings on the SEC's EDGAR system (public domain), fiscal years 2015–2016; `Data/corpus_item7/item7_2016_2017.qs2` is a `quanteda` corpus of 14,817 filings filed 2016–17 with document variables `cik`, `fyear`, `sic`. `prep_mdna.R` performs the paper's filters (13,309 fiscal-2015/16 filings; 13,120 after removing 189 duplicate firm-period observations; 11,491 after the 200-token length floor removes another 1,629, split 8,025/3,466). Simulation data are generated by seeded code — no external data.

**Deposit plan (JRSS-B policy).** The two data objects marked (DEPOSIT) are deposited on Zenodo with their own DOI. The raw corpus is excluded from the GitHub repository via `.gitignore` (file-size limit); the analysis-ready object is shipped here as well. The code repository is archived on Zenodo through the GitHub integration, one record per release, and the concept DOI is the one cited in the Data Availability Statement.

**Draft Data Availability Statement** (for the manuscript, before Acknowledgements):

> The Management's Discussion and Analysis texts analysed in Section 6 are extracted from Form 10-K filings publicly available on the SEC's EDGAR system (https://www.sec.gov/edgar). The extracted Item-7 corpus and the derived analysis-ready objects are deposited at [REPOSITORY], DOI [10.xxxx/xxxxx]. Simulation data are generated by seeded code. Complete replication code reproducing every exhibit in the paper and the Supplementary Material is available at https://github.com/contefranz/OpTop-R2-replication (archived at DOI [10.xxxx/xxxxx]); the methods are implemented in the R package OpTop (https://github.com/contefranz/OpTop).

## 8. Notes for referees

- `Code/README.md` is the pipeline manual: CLI grammar, every override key, engine details (WarpLDA vs VEM), monitoring of long runs, and per-experiment design notes with expectations. Where its inline examples use exploratory profiles, the commands in **this** README are the authoritative ones for the published exhibits.
- The null-discrepancy floor of Remark 1 (`min_null`, OpTop ≥ 0.14.1) is applied uniformly; excluded shares are reported in the outputs (`null_excl_share`, and the MD&A boilerplate index).
- Last verified inside this package on 2 September 2026 under OpTop 0.20.0 / NLPstudio 1.2.0 (R 4.6.1, macOS 27.0): `Rscript Code/tests_unit.R` passes U1–U16, and `./reproduce.sh exhibits` rebuilt every shipped csv/tex table byte-for-byte, both workbooks with identical sheet content, and all 23 shipped figure PDFs pixel-identically.
- Licensing: code under MIT; SEC filing text is in the public domain.
