---
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    theme: readable
    highlight: tango
---
# Simulation Pipeline

> **Read `../README.md` first.** This is the development manual of the pipeline (CLI grammar, override keys,
> engines, monitoring, design notes), written while the experiments were designed. The commands behind the
> manuscript are those of `../README.md` and `../reproduce.sh`. In particular, the published simulation design
> is K\* = 40, W = 5,000, J_train = 1,000, grid K = 10,...,100, WarpLDA with fitted priors alpha = 0.1 and
> beta = 0.01, ten training seeds (overrides `K_true=40 W=5000 K_grid=10:100:10 fit_method=WarpLDA alpha=0.1
> beta=0.01`); statements below about K\* = 10, VEM, 50 seeds or the `full` profile's defaults describe
> exploratory configurations, not the published runs. Run times for the published design are in `../README.md`.

Replication code for the simulation experiments of *"Goodness-of-Fit Indices
and Diagnostics for Topic Models"* (Grossetti & Lewis). Everything below runs
from the **project root**. Command grammar (all entry points):

```sh
Rscript Code/<script>.R <profile> [variant] [workers] [key=value ...] [label=name]
```

* `<profile>` — `smoke` (minutes, tiny corpus), `pilot` (small-scale but real),
  `full` (paper scale, $\sim 44$ h wall on 11 workers);
* a bare integer sets the **worker count** (else `OPTOP_WORKERS`, else cores $-$ 1);
* `key=value` pairs **override configuration parameters** (see below);
* `label=name` names the run; outputs go to `Results/Figures/<label>/…` and
  `Results/xlsx/Section5_tables_<label>.xlsx`. Without overrides the label is
  the profile; with overrides it auto-extends (e.g. `pilot_Kstar20_W5000`).

```sh
Rscript Code/make_all.R smoke                          # end-to-end exercise
Rscript Code/make_all.R pilot 8                        # pilot on 8 workers
Rscript Code/make_all.R full 11                        # paper scale
Rscript Code/make_all.R pilot 8 K_true=20 K_grid=5:40  # custom configuration
Rscript Code/run_E1.R full E1b 8                       # one robustness variant
Rscript Code/run_E1.R pilot 8 W=5000 L=500 label=smallvocab
Rscript Code/make_figures.R pilot W=5000 L=500 label=smallvocab  # same overrides!
```

### Override keys

Any scalar/vector field of the target experiment's config can be overridden;
unknown keys error with the valid list for that experiment. Grid/set-valued
keys accept three syntaxes:

```sh
K_grid=2:30                    # contiguous range
K_grid=5:100:5                 # stepped range (from:to:by)
K_grid=10,20,40,50,70,100      # arbitrary topic set
```

Adjacent gains and the $\varepsilon$-rule always compare *adjacent grid points* (Def. 1:
$K$ against $\mathrm{succ}(K)$), so sparse grids are handled correctly — the gain at
$K = 40$ on the set above is the improvement from 40 to 50 topics.

**Data-generating process** (all experiments)

| Key | Meaning | Baseline (pilot/full) |
|---|---|---|
| `K_true` | Number of topics in the generating LDA ($K^*$) | 10 |
| `W` | Vocabulary size | 10,000 |
| `J_train` | Training documents (models are fitted on these) | 1,000 |
| `J_eval` | Held-out evaluation documents (fresh draws, same DGP) | 500 |
| `alpha_DGP` | Dirichlet concentration of generated document–topic weights $\theta$ | 0.5 |
| `beta_DGP` | Dirichlet concentration of generated topic–word distributions $\varphi$ | 0.01 |
| `L` | Shorthand: fixed document length in tokens (sets `length_spec`). Does **not** affect E3, whose scenarios define their own length distributions | 1,000 |
| `seed_base` | Base of the deterministic seed schedule (change for an independent replication of everything) | 1970000 |

**Estimation and evaluation** (all experiments)

| Key | Meaning | Baseline |
|---|---|---|
| `K_grid` | Candidate topic numbers fitted on each training corpus; also defines the harmonized rare-word support (§3.2), so it is part of the index definition | 2:15 pilot / 2:20 full |
| `fit_method` | Estimator: `VEM` \| `Gibbs` (topicmodels) \| `WarpLDA` (text2vec) — see "Engines" below | `VEM` |
| `n_starts` | Restarts per fit; best training log-likelihood (topicmodels) / lowest training perplexity (WarpLDA) wins | 2 pilot / 3 full |
| `alpha` | Optional **fitted-model** document–topic prior; omit to retain the engine default. VEM fixes this value; Gibbs and WarpLDA also accept it. | engine default |
| `beta` | Optional **fitted-model** topic–word prior; omit to retain the engine default. Gibbs and WarpLDA accept it; VEM fails fast because topicmodels does not expose a fixed beta. | engine default |
| `c` | Harmonized-partition threshold ($\tau_j = c/L_j$; §3.2) | 1 |
| `metrics` | Discrepancy families to compute — `dev`, `chisq`, `se` (CLI can set a single one, e.g. `metrics=dev`; edit configs.R for subsets) | all three (E1/E3); `dev` (E2/E4) |
| `completion_prop` | Share of each evaluation document's tokens folded in under the completion target (rest is scored) | 0.5 |
| `eps_grid` | Tolerances $\varepsilon$ of the Definition-1 selection rule | 0.01, 0.005 |
| `sel_alpha` | One-sided level $\alpha$ of the $\varepsilon$-rule's upper bound and `optimal_topic()` | 0.05 |

**Replication counts** (which experiment uses which)

| Key | Meaning | Where |
|---|---|---|
| `S` | Corpus seeds (independent simulated corpora) | E1 (3/50), E3 (2/15 per scenario) |
| `S_train` | Training corpora whose fits are held fixed for conditional Monte Carlo | E2, E4 (2/10) |
| `R_eval` | Fresh evaluation sets per (training fit $\times J_{\mathrm{ev}}$) | E2 (60/300) |
| `J_ev_grid` | Evaluation-sample sizes at which coverage/size are measured | E2 (100,400 / 100,250,500) |
| `J_truth` | Size of the mega evaluation set that defines the conditional truth | E2 (4k/20k) |
| `R_null`, `R_power` | Evaluation replications for moment-test size / power | E4 (60,40 / 500,200) |
| `J_ev` | Evaluation-set size per moment-test replication | E4 (500) |
| `J_center` | Evaluation set used to estimate the conditional-truth centers of the moment tests | E4 (2k/20k) |
| `K_test` | Topic counts at which the moment tests run (auto-clipped to `K_grid`, $K^*$ always added) | E4 (5,10,15) |

**Experiment-specific tuning**

| Key | Meaning | Where |
|---|---|---|
| `gibbs_S` | **Deprecated, no effect** (the Gibbs robustness arm was removed); kept so old `gibbs_S=0` commands don't error | E1 |
| `B_strata`, `S_strata` | Frequency strata (Test 2) and fit strata (Test 3) counts | E4 (5, 5) |
| `min_docfreq` | Word filter for Test-3 strata and word-level tables (§3.8, Minor 5) | E4 (5) |
| `n_stop` | Size of the planted stopword pool in the contamination alternatives | E4 (50) |
| `c_grid` | Thresholds compared in the c-sensitivity arm | E5 (1, 5) |

`experiment` and `profile` are set by the positional arguments — don't
override them. List-valued fields (`length_spec`, E3's `scenarios`, E4's
`alternatives`, E5's `grids`) can't be expressed as `key=value`; edit
`Code/config/configs.R` for those.

**Shared overrides across `make_all.R`.** One override set is forwarded to every
stage, and experiments have different knobs — so a key is *applied where the
experiment uses it and ignored (with a one-line notice) elsewhere*. For example
`S=3` sizes E1/E3/E6 (`S`) and is skipped by E2/E4. **E2 and E4 size their replication differently:** use
`S_train=` (training corpora) and `R_eval=`/`R_null=`/`R_power=` (evaluation
replications), not `S`. A key valid for *no* experiment (a genuine typo like
`Ktrue=`) fails fast at startup with a nearest-match suggestion. When invoking a
single stage directly, its unknown keys are still ignored with a notice (the
driver reports which keys it used via the printed `tag=`).

`make_all.R` forwards overrides and label to every stage; when running stages
by hand, pass the **same** overrides to `make_figures.R`/`make_tables.R` so
they find the matching results.

Every stage is idempotent: fitted models are cached per (corpus, K) in
`Data/FITS/`, keyed by the (possibly overridden) DGP parameters, so
interrupted runs resume, the full profile reuses pilot fits, and different
configurations never collide.

### Monitoring long runs

Each driver first runs a **fit pool**: every (corpus $\times K$) fitting job across
all replicates is queued once, longest K first, so workers never
idle at replicate boundaries; the scoring pass afterwards reads the cache.
On entry the pool prints its composition — read it before walking away
(e.g. `WarpLDA 750 (K 10-80, 3 starts)`; a `VEM`/`Gibbs` line there means slow
fits that do NOT get WarpLDA's speed). A progress bar ticks per completed fit,
and each worker appends a line to the cumulative log:

```sh
tail -f Data/FITS/prefit.log      # [09:41:12] done K=80 WarpLDA (3 starts) in 214.3s
```

"Prefit" means *fitted before the scoring pass* — every logged fit is a
**genuine model** (the best of its restarts), cached to `Data/FITS/` and loaded
verbatim by the scoring pass; nothing is refit later, and every number in the
tables and figures traces back to exactly these fits.

**The scoring pass has its own log.** After the fit pool completes, each driver
enters a scoring pass that can run for a long time with a quiet console
(workers only relay output when a whole replicate resolves). Real-time
progress — replicate starts, per-stage ticks, and timings — streams to:

```sh
tail -f Data/scoring.log          # [15:12:40] E1 replicate 3/20: in-sample scored (412s)
```

Tip for custom big runs: at
$J \approx 10^4$ documents $\times$ $W \approx 2\times10^4$ terms the scoring-stage partition mask is
$\sim 0.8$ GB per corpus, so keep an eye on memory when many scoring workers run
in parallel.

**Regenerating a subset of figures/tables.** `make_figures.R` and
`make_tables.R` render every experiment by default; pass `exp=` to restrict:

```sh
Rscript Code/make_figures.R pilot exp=E3            # only E3 → Results/Figures/pilot/E3/
Rscript Code/make_figures.R pilot exp=E1,E3,E6      # a subset
Rscript Code/make_tables.R  pilot exp=E1
```

`exp=` selects *which* experiments; `label=` still controls *where* (the
`Results/Figures/<label>/` tree). E1 and E3 descriptive figures/tables show all
three discrepancy families (Deviance primary; Pearson and Squared-Error for
comparison); the Prop.-1(iii) gap decomposition, the $\varepsilon$-rule/selection, and the
§4 moment tests remain Deviance-based by construction.

**Discovering cached runs.** Data is keyed by the parameter tag, never by
`label=`, so re-plotting a finished run means reproducing its tag fields.
`list_runs.R` does that for you — it scans the caches and prints, per run, the
exact paste-ready command:

```sh
Rscript Code/list_runs.R                 # every cached run + its regen command
Rscript Code/list_runs.R E1 warplda      # filter: case-insensitive tag substrings
```

Each line reads a `Data/<exp>/<exp>_results_<tag>.qs2` cache and emits the
minimal `make_figures.R <profile> exp=<E*> …overrides…` that rebuilds its tag
(self-checked). Append `label=NAME` to choose the output folder; swap
`make_figures` $\to$ `make_tables` for that run's tables. Remember: the tag fixes
which *data* is read (`K_true`, `J_train`, `W`, `alpha_DGP`, `beta_DGP`,
`K_grid`, `fit_method`, and any explicit fitted `alpha`/`beta`), whereas
`label`, `S`, `S_train`, `gibbs_S`, and `J_eval` do
**not** enter the tag — so they never change which cache a command loads.

## Experiments

The paper's Section 5 presents **four studies**; the pipeline keeps one driver
per experiment (the merge is presentational):

| Paper study | Drivers | Main outputs |
|---|---|---|
| I — Held-out fit, selection & inference | `run_E1.R` + `run_E2.R` | F1/F1b/F2/F3; T1, T2a–e |
| II — Aggregation heterogeneity: documents & words | `run_E3.R` + `run_E6.R` | F6/F6b/F7/F11/F12; T6a/b, T9 |
| III — Moment-based specification tests | `run_E4.R` | F4/F5/F5b; T3–T5b |
| IV — Design sensitivity (appendix) | `run_E5.R` | F9; T7a/b |

Shared baseline DGP (paper §5.1): $K^* = 10$ true topics, $W = 10{,}000$ terms,
uniform document length $L = 1{,}000$, Dirichlet hyper-parameters $\alpha_{\mathrm{DGP}} = 0.5$,
$\beta_{\mathrm{DGP}} = 0.01$; $J = 1{,}000$ training and 500 evaluation documents (fresh draws from
the same DGP); LDA fitted by `topicmodels` VEM with multi-start selection;
harmonized-support threshold $c = 1$. Grid $K = 2{:}20$ (full) / $2{:}15$ (pilot).

### E1 — Held-out fit and topic-number choice (headline)

**What it does.** For each simulated corpus: fits the $K$ grid on the training
documents and scores three evaluation modes — in-sample (descriptive),
held-out-document *reconstruction* (fold in $\hat\theta$ from all of an evaluation
document's tokens), and held-out-token *completion* (fold in a binomial half,
score the other half) — computing all three $R^2$ families at Micro and Macro
level (§3.3–3.6), Prop-2 CIs, paired adjacent gains, and the Definition-1
selection $\hat K$ for $\varepsilon \in \{0.01, 0.005\}$. Comparators on the same corpora:
held-out perplexity and NPMI coherence (NLPstudio) and OpTop's original
`optimal_topic()` $\chi^2$ selector.
Replications: 3 seeds (pilot) / 50 seeds (full).
*(The former Gibbs and multi-start dispersion robustness arms were removed in
the 4-study consolidation; MC2/Minor-6 are addressed in the response letter.)*

**Expectation (theory $\to$ output).** Held-out curves rise steeply for $K < K^*$
and plateau after it, with completion strictly below reconstruction (a
stricter predictive target, §3.7); adjacent gains collapse to $\approx 0$ at $K^*$, so
$\hat K = K^*$. *Pilot: $\hat K = 10$ in 36/36 rule-combinations; perplexity and NPMI
also select 10; `optimal_topic()` rejects adequacy at every $K$ ($p \approx 0$,
df $\sim 10^6$) and its min-statistic fallback overshoots — the paper's
"binary tests vs effect sizes" motivation on display.*

**Outputs.** Figures F1/F1b (fit curves), F2 (gains + $\varepsilon$-rule),
F10 (variants, manual runs); tables T1 ($\hat K$ selection distribution incl.
comparators), T8 (runtime).

### E1b / E1c — Robustness variants (manual only)

Not part of `make_all.R`'s orchestration (no draft float, no referee item);
run by hand when needed: `Rscript Code/run_E1.R full E1b <workers>`.
E1b: $K^* = 20$, grid 5:40 (todo.txt item — is the elbow still sharp at higher
complexity?). E1c: $\alpha_{\mathrm{DGP}} = 0.1$ (near single-topic documents). Same scoring as E1.
**Expectation:** same qualitative pattern — held-out
elbow at $K^*$, $\hat K = K^*$ in most replicates, with noisier adjacent gains at
$K^* = 20$ (the paper's argument for a tolerance rule rather than exact visual
recovery). Output: F10 curves.

### E2 — Inference validation by conditional Monte Carlo

**What it does.** The paper's inference (Prop 2) is *conditional on the
training fit*, so each training grid is fitted once and only evaluation sets
are replicated: $R$ fresh evaluation samples per $J_{\mathrm{ev}} \in \{100, 250, 500\}$
(pilot: $R = 60$, $J_{\mathrm{ev}} \in \{100, 400\}$; full: $R = 300$), against conditional truths
($\mu$, true adjacent gains, true gap) estimated from one very large evaluation
set. Checks: empirical coverage of the 95% Macro CI, size of the paired
adjacent-gain test centered at the conditional truth, power against zero gain
for $K < K^*$, coverage of the delta-method gap CI (Remark 6), the $\hat K$ selection
distribution by $J_{\mathrm{ev}}$, and normal QQ plots of the $t$-statistics.

**Expectation.** Coverage $\approx 0.95$ at every $K$ and $J_{\mathrm{ev}}$; centered gain-test size
$\approx 0.05$; power $\to 1$ for $K < K^*$; true adjacent gain $\approx 0$ beyond $K^*$; gap CI
covers. This is referee MC1's non-negotiable ("coverage of the CIs in
Proposition 2, selection frequencies of $\hat K$"). *Pilot: coverage 0.88–0.97,
size 0.03–0.09, power $\approx 1$, gap coverage 0.91–0.96, $\hat K = 10$ in 480/480; the
measured true gain at $K^*$ was $-2\times10^{-5}$.*

**Outputs.** Figure F3 (coverage + QQ); tables T2a–T2e.

### E3 — Micro–Macro fit heterogeneity (upgraded old Experiment 2)

**What it does.** Three scenarios share the topic DGP and differ only in the
document population: **A** heterogeneous lengths (80% $\mathrm{Pois}(500)$, 20%
$\mathrm{Pois}(5000)$); **B** homogeneous lengths ($\mathrm{Pois}(1000)$) — the null case; **C**
equal lengths but two document groups with Dirichlet concentration $\alpha \in
\{0.2, 2\}$ — an *atypicality-driven* gap with no length heterogeneity at all.
Reports the held-out Micro–Macro gap with its delta-method CI (Remark 6), the
exact three-channel decomposition of Prop 1(iii) (length / atypicality /
interaction), subgroup Macro curves, document-level scatters (fit vs $L_j$ and
vs $\kappa_j$), and corrected length summary statistics (referee Minor 1).
Replications: 2 (pilot) / 15 (full) seeds per scenario.

**Expectation.** The gap is a *fit-heterogeneity* diagnostic, not a length
statistic (referee MC3): large and partly length-channel in A; small but
positive in B (pure atypicality — the honest null behavior, now with a CI);
large in C with the length channel **exactly zero** by construction. The
decomposition must sum to the gap at machine precision (identity). *Pilot:
gap at $K^*$ = 0.141 (A) / 0.036 (B) / 0.124 (C); C's length channel = 0.0000;
B's short/long subgroup Macro match to 3 decimals; scenario-A length SD
correctly $\approx 1800$.*

**Outputs.** Figures F6 (gap + CI), F6b (channel decomposition), F7
(scatters); tables T6a/T6b.

### E4 — Moment-based specification tests (§4): size, power, word diagnostics

**What it does.** Builds the §4 instruments from the training sample only
(Test 1 frequency contrast, $q = 1$; Test 2 five frequency strata, $q = 4$;
Test 3 five strata of the training word-level deviance $R^2$, $q = 4$) and runs
the Wald battery on replicated evaluation sets. **Size** (correct DGP):
rejection rates both *raw* ($H_0$: $\mu = 0$) and *centered at the conditional
truth* estimated from a large evaluation set — the raw null is rejected at
scale by any residual imbalance of the fitted model, e.g. VEM smoothing bias
(the paper's Remark 8), so the centered version isolates the calibration of
the Wald machinery, and effect sizes $\bar g$ (probability-mass units) are reported
throughout. **Power** (raw), by alternative $\times$ strength: doc-varying stopword
contamination (each document mixes its own random stopword subset — a shared
stopword distribution is just one extra LDA topic and gets absorbed; pilot
finding), eval-only contamination (clean training, contaminated evaluation),
burstiness (Dirichlet-multinomial), vocabulary drift (evaluation-side $\Phi$
perturbation), correlated topics (logistic-normal $\theta$). Companions: held-out $R^2$
across strengths, per-stratum $t$'s for one rejected case, and word-level
identification of the planted words by two criteria (mean held-out residual
$\bar e_w$; word-level $R^2$ with the precise filter of Minor 5: document frequency $\ge 5$
and baseline expected count $\ge 5$).

Alternatives (§5.2's five named cases, all implemented): doc-varying
contamination, eval-only contamination, burstiness, vocabulary drift,
correlated topics (logistic-normal $\theta$), and **document-group-specific
vocabulary** (`group_vocab`: $G$ document groups each over-weighting their own
vocabulary block — the natural showcase for the Test-3 fit-stratified and
Test-2 frequency-stratified instruments).

**Expectation.** Centered size $\approx 5\%$; raw size $\gg 5\%$ with tiny $\bar g$ ($\sim 10^{-7}$) —
report jointly per Remark 8. Power rises with strength; Test 3 catches
vocabulary-localized misspecification the frequency tests miss; the eval-only
contamination case is the headline: **overall $R^2$ barely moves while all three
tests reject with $\bar g \approx 100\times$ the null bias** — the "beyond overall fit"
message. Planted words top the residual ranking. *Pilot: centered size
T1 0.025–0.033, T2 0.033–0.067 (T3 0.12–0.17, explained by center-estimation
noise at $J_{\mathrm{center}} = 2000$; the full profile uses 20,000); eval-only
contamination at $w = 0.10$: $R^2$ 0.70 $\to$ 0.64 while all tests reject at 1.00;
planted words 50/50 in the top residual ranks; `group_vocab` rejects at power 1.0
under Tests 1–2.*

**Outputs.** Figures F4 (power), F5 (fit vs tests), F5b (word ranks); tables
T3 (size), T4/T4b (power, $R^2$), T5 (planted words), T5b (strata illustration).

### E5 — Design sensitivity (appendix)

**What it does.** One baseline corpus. (a) Grid extension: the harmonized
rare-word support is a union over the estimation grid (§3.2), so the SAME
fitted model's reported index changes when the grid grows — identical fits are
scored under nested grids (2:20 vs 2:50 in full) and the drift of $R^2(K)$ is
tabulated with min-bin shares (referee MC7). (b) Threshold sensitivity:
Deviance and Pearson indices at $c = 1$ vs $c = 5$ (referee MC8 — either $c = 5$
rehabilitates Pearson or it is demoted to an appendix).

**Expectation.** Grid-extension drift small on well-behaved corpora (the
indices are comparable across grids only up to the documented union effect —
the paper must state the (grid, $c$) convention); $c = 5$ curves sit higher
(smaller active support) and visibly steppier — evidence for the $c = 1$
default when deviance is primary. *Pilot: 2:15 $\to$ 2:22 extension moved the
reported indices imperceptibly; $c = 5$ higher and noisier for both families.*

**Outputs.** Figure F9; tables T7a/T7b.

### E6 — Word-level dual perspective (§3.8)

**What it does.** The frequency-space twin of E3: demonstrates the word-level
w-Micro vs w-Macro divergence (§3.8.2–3.8.3) and verifies **Lemma 2** (eq. 51,
the doc/word Micro-numerator commutativity on the unbinned support). Word-level
deviance indices are scanned over the whole $K$ grid (in-sample + held-out) under
the precise §3.8 filter (doc-freq $\ge 5$ & baseline expected $\ge 5$). Two scenarios:
**W1** correctly specified; **W2** high-frequency stopword contamination (the
example §3.8.3 names). Companion: per-word $R^2$ at $K^*$ vs document frequency with
the planted stopwords flagged. Replications: 2 (pilot) / 15 (full) per scenario.

**Expectation.** Both scenarios show w-Micro $>$ w-Macro (the model fits common
words better than the estimation-hard rare vocabulary — the frequency-space
length-bias analogue). W2's non-topic stopwords are mispredicted ($r^2 \approx 0$) yet
carry large baseline-discrepancy weight, so they **pull w-Micro down toward
w-Macro** (the "common words also mispredicted" regime) and populate the
worst-fit tail. Lemma 2 must hold at machine precision. *Verified ($W=3000$,
$K^*=8$): w-Micro 0.657 (W1) $\to$ 0.602 (W2) with w-Macro $\approx$ unchanged; planted
stopwords $5.7\times$ enriched in the worst-fit tail; Lemma-2 residual $8\times10^{-10}$.*

**Outputs.** Figures F11 (w-Micro/w-Macro over K, twin of F6), F12 (per-word
fit vs frequency); table T9 (gap at $K^*$ + Lemma-2 residual + planted-word
precision).

## Real-corpus application (MD&A, 10-K item 7)

Three scripts, separate from the E1–E6 orchestration, reuse the same modules
on the prepared corpus `Data/corpus_item7/item7_2016_2017.qs2` (quanteda
corpus, 14,817 filings filed 2016–17; docvars incl. `cik`, `fyear`, `sic`).
The two **fiscal years** (the year an MD&A discusses) are **pooled** into one
corpus with a stratified random train/held-out split:

```sh
Rscript Code/prep_mdna.R  y1=2015 y2=2016                  # corpus -> dtms (~minutes)
Rscript Code/run_mdna.R   y1=2015 y2=2016 workers=6        # fits + analyses A-D
Rscript Code/make_mdna_outputs.R y1=2015 y2=2016           # M1-M4 figures + tables
```

* `prep_mdna.R` — fiscal-year filter, (cik × period) dedupe, parallel
  tokenization via `NLPstudio::tokenize_corpus()` (`workers=4`;
  punct/symbols/numbers out, ≥ 4 chars, lowercase, English stopwords out),
  train/eval split of y1 (`train_frac=0.7`, seeded), training-vocabulary
  pruning (`w_mindoc=5`), y2 projected on the training vocabulary (OOV mass
  reported), length floor (`len_floor=200`), FF12 industry from SIC. Saves
  `Data/MDNA/mdna_prep_<y1>_<y2>.qs2`.
* `run_mdna.R` — WarpLDA (3 starts, **engine-default priors**) on a coarse
  grid (`K_grid=25:200:25`) with automatic refinement around the elbow;
  then (A) fit/selection incl. comparators, (B) consistency battery at K̂
  (Micro/Macro + Prop-2 CIs per family, gap + Prop-1(iii) channels,
  **boilerplate index** = null-floor exclusion share, also by FF12),
  (C) word-level dual + worst-fit vocabulary, (D) §4 moment tests **on the
  held-out set** with training instruments (raw + effect sizes; no
  conditional truth on real data) + signed residual ranking (over-observed
  vocabulary).
* `make_mdna_outputs.R` — the paper set, under
  `Results/Figures/mdna_<y1>_<y2>/MDNA/` +
  `Results/xlsx/MDNA_tables_<y1>_<y2>.xlsx`. Figures: **R1** (fit curves,
  3 protocols), **R1b** (Micro–Macro gap + delta-method CI), **R1c**
  (adjacent gains + ε rule), **R4** (moment diagnostics: per-stratum mean
  moments at K̂, and effect sizes over the tested grid); tables **R1**
  (selection), **R2** (battery at K̂), **R3** (vocabulary); appendix A1
  (word dual), A3 (fit vs length), A2 (boilerplate by industry).
  `families=dev` (default) renders Deviance-only figures — family agreement
  is certified by the R1 selection table; `families=all` restores the
  three-family facets for an appendix variant. For R4's over-K panel run the
  driver with `tests_K=all` (moment tests at every grid point, roughly one
  extra scoring pass; default `tests_K=hat` tests K̂ ± 1 only).

**Smoke path** (run this first): pass `sample_n=200` to all three scripts,
plus a small grid to the driver (`K_grid=10:30:10 refine_span=0`) — 200 docs
per fiscal year end-to-end in minutes. Monitoring: `tail -f Data/scoring.log`
as usual.

`make_mdna_outputs.R` also takes `in_suffix=` (read a re-scored object, e.g.
`in_suffix=_rev1`), `out_label=` (explicit output label) and `K_ref=` (the
reference fit whose diagnostics are drawn; must equal the object's `K_hat`).
Figures, csv, tex **and the workbook** follow the label, so a run on a re-scored
object never overwrites the baseline exhibits.

## September 2026 revision tooling (evaluation only — no model is estimated)

Two defects were fixed and the selection rule was revised; every recomputation
reuses the cached fits in `Data/FITS`.

* **Test-3 instruments.** The fit-stratified instruments were built from the
  strata of the *smallest* tested K for every K (a `data.table` scoping defect:
  a closure argument named like a column is shadowed by the column —
  `word_tr[K == K & …]` in `run_mdna.R`, `word_tr_l[list(K), ws]` in
  `run_E4.R`). Both drivers now call `make_instruments_by_K()`
  (`R/utils_moment_tests.R`); gate **U17** tests the drivers' code path. Tests
  1–2, every fit index and Test 3 at the smallest tested K were never affected.
* **Selection rules** (`R/utils_inference.R`): `paired_gains_all()`,
  `select_k_total_gain()` (primary: Bonferroni over all m(m−1)/2 pairs),
  `select_k_adjacent()` / `select_k_adjacent_simultaneous()` (local companion),
  `select_k_all_rules()` (all three, incl. the original pointwise rule, now
  exploratory). "None certified" is `NA`, never the grid maximum; optional
  cluster-robust standard errors (`cluster=`). Gate **U18**.
* **Held-out word-level Deviance** uses the Poisson-form baseline deviance
  `2 Σ_j [N log(N/B) − (N − B)]` on the *scored* tokens
  (`.word_null_dev_poisson`, gates **U19**, **U25**). OpTop ≤ 0.20.1 omits the
  linear term in the word-level null (zero in-sample, non-zero held-out);
  0.20.1 raises a classed warning (`optop_word_null_baseline`) and leaves the
  kernel unchanged, so the pipeline keeps its own null. The warning is muffled
  only in the branch that replaces that null (gate **U20**).
* **One convention, applied once** (`R/utils_revision.R`). Every new result
  object is stamped (`run_meta$scoring`: `word_null_convention =
  "poisson_scored_tokens"`, OpTop version and commit). Post-processing converts
  a word-level null only when the input is a *registered* pre-revision object,
  identified by content hash (`LEGACY_WORD_NULL_INPUTS`); a stamped object
  passes through unchanged and anything else is an error (gate **U23**).
* **c and δ are separate.** `run_mdna.R` takes `c_part=` (support threshold)
  and `min_null=` (discrepancy floor); design comparisons over c hold δ = 1 and
  use `refine_span=0`, so the same candidates enter every harmonised support
  (gate **U21**).
* **Intermediates kept.** MD&A and E6: per-word fitted and null discrepancies at
  every K. E4: conditional centres and, per replication, the projected-moment
  summaries (n, mean, covariance, centre) from which every Wald statistic can be
  recomputed. E1/E2: all-pairs gains and the selections under the three rules.
  E1/MD&A/E5: support-resolution rows (`support_resolution()`, gate **U24**).
* **Checkpoints** (`revision_checkpoint()`): with an `out_suffix`, drivers store
  each completed seed / battery / evaluation block under
  `Data/Checkpoints/<suffix>/` and reuse it only when configuration, package and
  *scoring* code are identical; a mismatch is an error (gate **U22**).
* **No-fit guard.** `options(optop.no_fit = TRUE)` or `OPTOP_NO_FIT=1` turns a
  fit-cache miss into an error (`prefit_pool()`), so a drifted configuration can
  never start estimating models silently.
* **Driver options.** `cfg_from=<results.qs2>` reuses the configuration stored in
  a cached result object (identical corpus signature, seeds and fit spec, hence
  identical fit-cache keys); `out_suffix=` namespaces every output. Nothing
  produced for the revision overwrites a pre-revision object: outputs carry
  `_rev1` (fast phase) or `_rev2` (batch).

```sh
Rscript Code/tests_unit.R                                  # U1-U25 (~10 s)
Rscript Code/postprocess_revision.R                        # cache-only, 19 gates (~1 min)
Rscript Code/rescore_revision.R stages=tests,restarts,lemma,wordfloor   # ~15 min
Rscript Code/make_revision_tex.R                           # Results/tex/rev1_*.tex
Rscript Code/make_revision_figures.R                       # Results/Figures/*_rev1/
# the long evaluation-only re-runs: resumable, verified, fit cache hashed
bash Code/run_revision_batch.sh 8 A                        # MD&A, E2, E1, E6, diagnostics (~15 h)
bash Code/run_revision_batch.sh 8 B                        # E4 both arms, E1 10 x 10 (~18 h)
bash Code/run_revision_batch.sh 8 final                    # post-processing + acceptance
bash Code/run_revision_batch.sh 8 all _rev2 --dry-run      # what is complete, what would run
```

**The batch** (`revision_batch.py`, launched by `run_revision_batch.sh`). A stage
is complete only if it exited 0, its declared outputs exist, and its identity —
command, input hashes, package, and the code that determines its result — is the
one on disk. Identity is per stage: the shared modules (`R/*.R`),
`config/configs.R` and the script the stage runs, so repairing an exhibit
script, or a different driver, never invalidates a night of scoring. Unit gates and the toy smoke run
(`revision_smoke.sh`: every driver, the resume path, `cfg_from=` under the no-fit
guard) are prerequisites and re-run automatically when the scoring code or the
package changes. Every production stage runs with `OPTOP_NO_FIT=1`; the sha256
of every cached fit is compared with `Results/csv/fit_manifest_rev1.csv` before
and after. Completed, verified, failed, blocked and skipped stages are reported
separately (`Results/revision_rev2/summary.json`, `stages.json`, one log per
stage); the runner never reports completion of stages it did not run.

| id | key | what |
|---|---|---|
| 2 | `mdna` | MD&A rescoring over K = 10–200: word curves and raw word discrepancies, tests at K ∈ {40,50,60,170,180,190}, support resolution |
| 3 | `mdna_c05`, `mdna_c2`, `mdna_g100` | c ∈ {0.5, 2} at δ = 1 on the same 20-model grid; grid truncated at 100 |
| 4 | `e2_base`, `e2_dgp2` | 20,000 reference documents per training fit, reference SEs, three rules |
| 5 | `e1_base`, `e1_dgp2` | ten seeds: all-pairs gains, three rules, support resolution |
| 6 | `e6` | word study with corrected nulls, raw discrepancies, independent Lemma S1 check |
| 7 | `e4_base`, `e4_prior` | K-specific Test-3 strata; centres and moment summaries kept |
| 8 | `selreps` | 10 × 10 evaluation corpora, three rules |
| 9–11 | `e5_resolution`, `unbinned`, `restarts` | `run_revision_diagnostics.R` (E5 support resolution; unbinned protocol-matched comparator with φ floored at 10⁻¹², 10⁻¹⁴, 10⁻¹⁰ and renormalised); restarts on the support of the production grid and the restarts |
| 12–13 | `post_*`, `acceptance` | `postprocess_e4.R`, `postprocess_e2_gap.R`; `finalize_revision.R` → `acceptance.json` |

* `postprocess_revision.R` — selections under the three rules (MD&A, E1, second
  DGP), firm-cluster standard errors, δ and design sensitivity, residual mass by
  frequency group, E2 coverage calibration against a noisy reference, planted
  words. Stops on a failed regression gate. → `Results/csv/*_rev1.csv`,
  `Data/MDNA/revision_postprocess_rev1.qs2`.
* `rescore_revision.R` — corrected MD&A moment tests at K ∈ {40, 50, 60, 170,
  180, 190} with iid and firm-clustered covariance, restarts on a common
  support, the Lemma-S1 cross-path check, word-level conventions. Its first gate
  forces the K = 10 strata and must reproduce the *published* Test 3 exactly
  before any corrected value is written. → `Data/MDNA/mdna_rescore_rev1.qs2`.
* `finalize_revision.R` — acceptance checks on the batch outputs, read from the
  explicit paths the runner recorded: stamps, configurations, the held-out
  baseline deviance of the word "condition" recomputed from the definition
  (9562.168059; legacy 8771.558890), E6 Lemma S1 across code paths, untouched
  scores against the published caches, corrected MD&A Test 3 against the
  verified values, fit hashes.

Provenance of every number in the manuscript: `../README.md`, Section 6 (exhibit map) and `Results/manuscript_files.csv`.

## Layout

```
Code/
  R/                  modules (sourced via R/source_all.R)
    utils_io.R        paths, run tags, config hashes, qs2 caches, seeds log
    utils_dgp.R       LDA DGP + misspecification injectors
    utils_fit.R       topicmodels VEM/Gibbs fitting, fold-in, token splits,
                      pseudo `nlp_topic_fit` objects (OpTop's official adapter)
    utils_heldout.R   §3.7-exact held-out partition/baseline + scoring
    utils_inference.R Prop-2 CIs, paired gains, ε-rule, gap SE, Prop-1iii
    utils_moment_tests.R  §4 instruments + Wald tests (raw and centered)
    utils_comparators.R   perplexity/NPMI (NLPstudio), optimal_topic (OpTop)
    utils_revision.R  scoring provenance, word-null convention guard,
                      checkpoints, support resolution
    theme_paper.R     shared ggplot style/encodings
  config/configs.R    all parameters; profiles smoke / pilot / full
  run_E1.R … run_E6.R experiment drivers (idempotent, cache-backed)
  run_E1_selreps.R    Study-I selection on S x R_sel replicated evaluation
                      corpora (evaluation-only; cached fits)
  run_mdna_restarts.R across-restart dispersion of the MD&A held-out index
  postprocess_e4.R    size-adjusted power + contamination-refit seed profile
                      from a saved E4 object (no fits)
  postprocess_e2_gap.R gap-CI coverage diagnosis from a saved E2 object
  run_tonight.sh      one-shot robustness batch (stages 0-7, resilient,
                      resumable: `bash Code/run_tonight.sh <workers> <from>`)
  postprocess_revision.R  Sept-2026 revision: cache-only recomputations (gated)
  rescore_revision.R      Sept-2026 revision: short re-scoring on cached fits
  make_revision_tex.R     Sept-2026 revision: table bodies -> Results/tex/rev1_*
  make_revision_figures.R Sept-2026 revision: changed figures -> Figures/*_rev1
  run_revision_batch.sh   Sept-2026 revision: launcher of the evaluation-only
                          batch (`bash Code/run_revision_batch.sh <workers> <A|B|final|all>`)
  revision_batch.py       stage table, identities, verified resume, fit hashing
  revision_smoke.sh       toy end-to-end run of every batch driver
  run_revision_diagnostics.R  E5 support resolution; unbinned comparator
  finalize_revision.R     acceptance checks on the batch outputs
  tests_unit.R        correctness gates (run first; non-zero exit on failure)
  make_figures.R      F1–F12 from saved results only
  make_tables.R       T1–T9 (tinytable tex + CSV + one Excel workbook)
  make_all.R          orchestrator
  list_runs.R         list cached runs + their exact re-plot commands
  99_LEGACY/          superseded scripts (working tree only; not shipped)
```

Outputs (everything under `Results/` and `Data/`):

* `Data/<E*>/….qs2` — tidy results with embedded config + package versions
  (filenames carry the full parameter tag, so configurations never collide);
* `Data/FITS/` — per-(corpus, K) fit cache shared across experiments/profiles;
* `Results/Figures/<run_label>/E1..E5/` — one self-contained figure tree per
  configuration (default labels: `smoke/`, `pilot/`, `full/`); `OLD/` and
  `NIPS/` hold relocated legacy material;
* `Results/xlsx/Section5_tables_<run_label>.xlsx` — **one workbook per
  configuration, one sheet per table** (browse this first);
* `Results/tex/`, `Results/csv/` — LaTeX fragments and CSV copies, suffixed
  by run label;
* `Results/csv/seeds.csv` — every seed used (replication package).

Note: the legacy NIPS real-corpus scripts (`run_nips_section3.R`,
`prepare_nips_corpus.R`, `utils_sim.R`) are out of this pipeline's scope and
are not part of the replication package; they live only in the project's
working tree and write to a top-level `Figures/NIPS` path.

## Engines

Fitting goes through the `NLPstudio::fit_topic_model()` facade; the returned
`nlp_topic_fit` is the model class OpTop's adapter supports natively, so the
whole scoring path is engine-agnostic. `fit_method=WarpLDA` switches to
text2vec's WarpLDA — measured **$\sim 100\times$ faster** than VEM on this machine
($K=20$, $J=500$, $W=5000$: 2 s vs minutes) — with two caveats:

* WarpLDA is collapsed MCMC with **fixed priors** by default
  (`doc_topic_prior = 50/K`, `topic_word_prior = 0.1`; no $\alpha$ estimation),
  although explicit fitted `alpha`/`beta` overrides are passed through to those
  two priors. Its $\hat\theta$ is from the final sampler state and can be noisier
  than VEM posterior means for very short documents;
* its backend object cannot be serialized, so WarpLDA fits are cached slim
  ($\theta/\varphi$ matrices only) and held-out fold-in for cached fits uses a fixed-$\varphi$ EM
  (validated against `topicmodels::posterior()` in `tests_unit.R` U8).

Recommendation: keep the **paper's headline numbers on VEM** (matches §5.1
and the Remark-8 narrative in E4); use WarpLDA for exploration, large-K
sweeps, and as an extra estimator-robustness overlay.

## Design notes

* **Held-out protocol (§3.7).** Training baseline $\hat\pi^{\mathrm{tr}}$ and training topics
  only; $\hat\theta$ for evaluation documents by fold-in (`topicmodels::posterior`);
  *reconstruction* folds in the full document, *completion* folds in a
  binomial half and scores the other half. The harmonized rare set uses
  $\min(\hat\pi^{\mathrm{tr}}(w),\ \min_K \hat\imath_{jw}) < c/L_j$ with $L$ from the scored tokens.
  Since OpTop 0.14.0 `optop_make_partition(pi_glob=)` accepts the external
  training baseline, so `make_heldout_partition()` delegates the partition to
  the package's compiled cores; the pure-R reference implementation is kept as
  `.make_heldout_partition_ref()` and native≡reference equality is enforced by
  `tests_unit.R` (U15), alongside the in-sample crosscheck (U1).
* **Null-discrepancy floor (OpTop 0.14.1).** A document whose harmonized
  support collapses into the pooled min bin at resolution $\tau_j = c/L_j$ has
  observed $\equiv$ baseline by construction, so $D_j(\mathrm{null}) \approx 0$ in **every**
  discrepancy family at once and $1 - D_j(K)/D_j(\mathrm{null})$ is unbounded —
  a handful of such documents destroy the Macro curves and their CIs (the
  completion protocol at low-atypicality design points is the exposed case;
  Micro, being the $D_{\mathrm{null}}$-weighted mean, is immune). Since OpTop 0.14.1 the
  index functions exclude documents with $D_j(\mathrm{null}) < c$ (`min_null`,
  default = the partition constant) uniformly across Dev/Pearson/SE; the
  pipeline reports the excluded share as `null_excl_share` next to `J_pos`,
  and the paper's Prop-2 conditioning event sharpens from
  $\{D_\mathrm{null} > 0\}$ to $\{D_\mathrm{null} \ge c\}$.
* **Why OpTop 0.14.0's held-out functions are not adopted wholesale.**
  `optop_index_holdout()` folds $\hat\theta$ in via *live* backend objects, but this
  pipeline caches WarpLDA fits *slim* (the C++ pointer does not serialize), so
  scoring must go through the local fixed-$\varphi$ EM fold-in; the package also has
  no completion token-split scorer, no word-level held-out indices, and
  `optop_moment_test()` has no centered/conditional-truth null (E4's Remark-8
  design). The local scorer is validated against `optop_index_holdout()` on
  live fits in U15 (agreement at $10^{-8}$).
* **No namespace patching.** Held-out $(\hat\theta, \hat\Phi)$ pairs ride OpTop's exported
  index engine as minimal objects of class `nlp_topic_fit` (an officially
  supported adapter class).
* **Conditional Monte Carlo (E2/E4).** Inference in the paper is conditional
  on the training fit, so coverage/size studies fix the fitted grid and
  replicate evaluation sets only — hundreds of replications without refits.
* **Fit cost** (reference machine, $J=1000$, $W=10000$, $L=1000$, VEM one start):
  $K=5$: $\approx 27$ s, $K=10$: $\approx 92$ s, $K=20$: $\approx 330$ s, $K=40$: $\approx 19$ min; scenario-A corpora
  $\approx 3\times$ slower. The prefit pass pools every (corpus $\times K$) job, longest first.

* **OpTop ≥ 0.15 compressed partition ("format 2").** `optop_make_partition()`
  returns per-document 0-based indices of the non-rare words instead of the
  dense `rare_mask`; `.densify_partition()` (`utils_heldout.R`) rebuilds the
  dense mask the pipeline consumes, a no-op for ≤ 0.14 objects. OpTop 0.16
  removed the deprecated `ztest=`/`reopt=` index arguments and 0.19 renamed
  `optimal_topic()` to `optop_select()` (delegating alias kept); the pipeline
  uses neither removed argument. Gates U1/U2/U15/U16 enforce agreement with
  the frozen pure-R reference at 1e-10 under both 0.14.1 and 0.20.0; OpTop's
  changelog documents two bit-level changes of at most ~1e-12 (0.15.0, 0.17.0).
  The 0.14.1 null-discrepancy floor never binds under reconstruction or
  in-sample on the paper's designs (E1 stamps: `null_excl_share` = 0 there,
  0.11–0.20 under completion), so the two caches produced under 0.14.0 (E4,
  E6; reconstruction/in-sample only) are unaffected.

## Dependencies

R $\ge$ 4.6 with: OpTop (0.20.1, pinned by `install.R`; 0.20.0 and 0.14.1 also verified for everything but gate U20),
NLPstudio (1.2.0, pinned; 1.1.1 also verified), topicmodels, quanteda,
data.table, Matrix, future.apply, ggplot2, patchwork, qs2, tinytable, writexl,
here, digest, MASS.
