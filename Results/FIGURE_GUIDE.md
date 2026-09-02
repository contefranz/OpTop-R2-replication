# Section 5 figures — an intuition guide

A plain-language tour of the simulation figures. Each entry has
three beats: **what you see**, **how to read it**, and **the takeaway**. The
E1–E5 figures are shown in `Results/Figures/pilot/` (VEM engine) and the
complete set including **E6** in `Results/Figures/warptest/` (the WarpLDA run,
which ran all six experiments). Numbers quoted are from these *small pilot* runs
(few documents, noisy on purpose); the full run sharpens them but the shapes are
already the story. The paper's primary index is **Deviance** — read it first in
every multi-panel figure; Pearson and Squared-Error are shown alongside for
contrast.

**Paper ↔ pipeline map.** Section 5 presents **four studies**; the pipeline's
experiment drivers (and the figure directories) keep their E-labels:

| Paper | Pipeline drivers | Figures |
|---|---|---|
| Study I — held-out fit, selection & inference | E1 + E2 | F1, F1b, F2, F3 |
| Study II — aggregation heterogeneity: documents & words | E3 + E6 | F6, F6b, F7, F11, F12 |
| Study III — moment-based specification tests | E4 | F4, F5, F5b |
| Study IV — design sensitivity (appendix) | E5 | F9 |

---

## How to read any of these figures

**The core question — `R²_D(K)`.** Imagine two reference models. The *floor* is
the "no-topics" baseline: every document is assumed to look like the
corpus-average bag of words. The *ceiling* is a perfect reconstruction of the
observed counts. A K-topic model sits somewhere between. `R²_D(K)` is **the
fraction of that floor-to-ceiling gap the model closes** — 0 means "no better
than assuming every document is average," 1 means "perfect." It's the topic-model
analogue of regression R². The **dashed vertical line** in most plots marks the
*true* number of topics used to generate the data, `K* = 10`.

**Three discrepancy families** (the "D" in `R²_D`) — three ways to measure how
far fitted counts are from observed:
- **Deviance** — likelihood-based, the one LDA effectively optimizes. Best-behaved and the paper's primary.
- **Pearson (χ²)** — connects to classical χ² adequacy tests, but its *level* is unstable when expected counts are tiny; read its shape, distrust its exact height.
- **Squared-Error** — raw squared count errors; dominated by long documents and frequent words. A geometric complement, not the headline.

**Micro vs Macro aggregation** — two ways to average document-level fit into one
number. **Micro** pools all documents (so long, unusual documents dominate);
**Macro** is the plain average over documents (every document counts equally).
When they disagree, the fit is *concentrated* somewhere — that gap is itself a
diagnostic (see E3).

**Three evaluation modes** — a ladder from optimistic to honest:
- **In-sample** (solid line, ● circle) — score the same documents the topics were trained on. Descriptive only; always the rosiest.
- **Held-out (reconstruction)** (long-dash, ▲ triangle) — new documents; infer each one's topic mix from *all* its words, then score it. The document "sees itself."
- **Held-out (completion)** (dotted, ■ square) — new documents; infer the topic mix from *half* the words and predict the *other* half. The model can't cheat — this is the genuine predictive curve, closest to held-out perplexity.

Completion always sits **below** reconstruction; the gap between them is the
"optimism" of letting a document inform its own score.

**Visual conventions.** Line style *and* point shape both encode the series
(so overlapping curves stay legible in black-and-white); facet columns are the
discrepancy family or the scenario; shaded ribbons are confidence intervals or
the spread across random seeds.

---

## E1 — Does the fit tell you how many topics?

### `F1_heldout_curves` — the headline
**What you see.** `R²_Macro` climbing with K, in three panels (Deviance,
Pearson, Squared-Error); within each panel, the three evaluation modes.

**How to read it.** Look for the **elbow**: steep gains up to `K* = 10`, then a
flat plateau. The dotted-square (completion) line runs a notch under the
long-dash-triangle (reconstruction) line everywhere. One caveat since OpTop
0.14.1: the completion Macro may exclude near-baseline documents via the
null-discrepancy floor (`D_null < c`) — the share sits in `e1_summary`'s
`null_excl_share` column and in table T1b, and it is the *same documents at
every K* (the null discrepancy doesn't depend on the fitted model), so it
shifts the level, never the elbow.

**Takeaway.** The elbow lands exactly at the true `K* = 10` in all three
families — additional topics past 10 buy almost nothing. The completion-below-
reconstruction gap is the honest "predicting unseen words is harder than
re-describing seen ones." That gap is *cleanest in Deviance*, artificially small
in Pearson, and *largest and noisiest in Squared-Error* (its wide ribbon on the
right) — a one-glance illustration of why Deviance is the primary index.

### `F1b_micro_curves` — the Micro companion
Same as F1 but pooled (Micro) instead of averaged (Macro). Here the two look
almost identical because E1's documents are all the same length — when nothing
makes some documents dominate the pool, Micro and Macro agree. (E3 is where they
part ways.)

### `F2_adjacent_gains` — the stopping rule in action
**What you see.** For each step up the K-ladder, the *extra* `R²` that one more
topic delivers, with a one-sided upper error bar; the panels are the two
held-out modes × three families. Horizontal dashed lines mark tolerances
ε = 0.01 and 0.005.

**How to read it.** Walk left to right; the first K whose entire error bar drops
**below** the ε line is the chosen `K̂`. The dot-dash vertical line marks it.

**Takeaway.** Gains are large and clearly above ε up to K = 10, then collapse
into the noise band — so the rule selects `K̂ = 10`, unanimously across ε, both
targets, and all families (36/36 in the pilot). This is the paper's point that
you don't eyeball the elbow — you test whether the *next* topic is worth it.

### `F8_estimator_robustness` — REMOVED
*(The Gibbs/multi-start robustness arms were removed from E1 in the 4-study
consolidation; the pipeline no longer produces this figure. Old copies under
`pilot/E1/` are historical. The MC2/Minor-6 point — VEM and Gibbs trace the
same curve, restart spread negligible — is handled in the referee response.)*

---

## E2 — Can you trust the confidence intervals?

### `F3_coverage_qq` — are the error bars honest?
**What you see.** Top: empirical coverage of the nominal 95% interval for the
average held-out fit, plotted against K, in panels for different evaluation-set
sizes `J_ev`. Bottom: normal QQ-plots of the underlying t-statistics.

**How to read it.** Coverage should hug the dashed 0.95 line; the QQ points
should lie on the diagonal (i.e., the statistic really is standard-normal).

**Takeaway.** Coverage sits at 0.88–0.97 (within Monte-Carlo error of 95% at
this pilot size) and the QQ points track the diagonal — the paper's central
limit theorem and its intervals are **calibrated**. This is the evidence that
the whole inferential apparatus (confidence intervals, the paired test behind
`F2`) can be trusted, which is the paper's core promise.

---

## E3 — Who is the "average fit" hiding?

Three scenarios share the same topics but differ in the documents: **A** mixes
short and very long documents; **B** makes them all the same length (the null);
**C** keeps lengths equal but makes some documents sharply single-topic and
others diffuse (different "typicality," no length difference).

### `F6_gap_ci` — the Micro–Macro gap
**What you see.** `R²_Micro − R²_Macro` over K, in a 3×3 grid (family rows ×
scenario columns), each with a delta-method 95% CI ribbon.

**How to read it.** A gap near zero means the pooled and averaged fit agree
(fit is spread evenly). A positive gap means the pooled index is being lifted by
a minority of documents that fit especially well.

**Takeaway.** Scenario A opens a large gap (≈0.14 in Deviance) — the long
documents, which dominate the pool, are fit better than the typical short one.
Scenario B stays small but *not exactly zero* (≈0.036) — and now it has a
confidence band, so we can say the residual gap is a real, quantifiable
"typicality" effect rather than hand-wave it as "close to zero." The gap is a
*fit-heterogeneity* signal, and — reading down the family rows — Squared-Error
shows the biggest Scenario-A gap of all (its quadratic sensitivity to length).

### `F6b_gap_decomposition` — *why* the gap exists
**What you see.** The Deviance gap split into three stacked contributions —
**length**, **atypicality**, and their **interaction** — with the total gap as
a line on top.

**How to read it.** Which colored band dominates tells you the *source* of the
heterogeneity.

**Takeaway.** In Scenario A the length band is large; in Scenario C the length
band is **exactly zero** (by construction the documents are all the same
length) and the whole gap is the atypicality channel. This is the decisive
picture that the Micro–Macro gap is **not merely a document-length artifact** —
it responds to how *unusual* a document is, not just how *long*.

### `F7_doc_scatters` — the mechanism at the document level
**What you see.** One dot per evaluation document: left, fit vs document length
(Scenario A); right, fit vs `κ` = how far the document's word distribution sits
from the corpus average (Scenario C). A binned median line runs through each.

**How to read it.** An upward trend means "this property predicts good fit."

**Takeaway.** Left: longer documents fit better (more tokens → the topic mix is
pinned down more precisely). Right: more atypical documents fit better even at
equal length. Together they *are* the two channels of `F6b`, made concrete at
the level of individual documents.

---

## E4 — Catching problems the R² can't see

The R² answers "how much fit," not "what's still wrong." E4 plants known defects
and shows the moment-based specification tests catching them.

### `F4_power_curves` — do the tests fire when the model is wrong?
**What you see.** Rejection rate of the three tests (frequency-contrast,
frequency-strata, fit-strata) as a misspecification is dialed up, one panel per
defect type.

**How to read it.** Curves should rise from the nominal 5% (the dotted "size"
reference) toward 1.0 as the defect strengthens; a curve stuck near 0.05 is a
test blind to that defect.

**Takeaway.** Power climbs with defect strength, and the *fit-strata* test in
particular catches vocabulary-localized problems the frequency tests miss —
demonstrating the tests do what a single R² number cannot: localize *where* the
model fails.

### `F5_fit_vs_tests` — the punchline
**What you see.** On the same x-axis (defect strength): the overall held-out
`R²` and the moment-test rejection rate.

**How to read it.** Watch the two lines diverge — R² barely dips while rejection
shoots to 1.

**Takeaway.** Under evaluation-side contamination the overall fit slips only
slightly (≈0.70 → 0.64) yet the tests reject essentially always. This is the
paper's "beyond overall fit" message in one image: **a model can look fine on
aggregate R² and still be detectably misspecified.**

### `F5b_word_ranks` — naming the culprits
**What you see.** Every vocabulary word ranked by its average held-out residual
(how much its observed frequency exceeds what the model predicts); the planted
"bad" words are highlighted.

**How to read it.** If the highlighted points cluster at the top of the ranking,
the diagnostic successfully fingered the injected words.

**Takeaway.** The planted contamination words sit right at the top of the
residual ranking (top ~50 in the pilot) — so the word-level diagnostic doesn't
just say "something's wrong," it hands you the specific vocabulary to inspect.

---

## E5 — How sensitive are the knobs?

### `F9_design_sensitivity`
**What you see.** Left: the Deviance index computed under two nested topic grids
(does extending the grid move the curve?). Right: the index at two values of the
rare-word threshold `c` (1 vs 5), for Deviance and Pearson.

**How to read it.** Overlapping curves = the choice doesn't matter; separated
curves = it does, and you must report it.

**Takeaway.** Grid extension barely moves the index (reassuring — the harmonized
support is stable). The threshold `c` matters more: `c = 5` runs higher and
choppier than `c = 1`, which is why the paper fixes `c = 1` for the
deviance-primary analysis and treats `c` as a reported sensitivity rather than a
silent default.

---

## E6 — Where in the vocabulary does the model fail?

E3 asked which **documents** the average hides; E6 asks the same question
word-by-word, on the correctly specified corpus (**W1**). *(The former W2
scenario — a shared stopword block injected on both sides — was dropped from
the paper: a shared stopword distribution is just one extra topic, so the
refit model absorbs it and the word-level indices barely move. That is the
same absorption finding that led E4 to use doc-varying/eval-side
contamination. The pipeline still computes W2; the figures show W1.)*

### `F11_word_micro_macro` — the frequency-space twin of the Micro–Macro gap
**What you see.** Word-level Deviance R² over K: **w-Micro** (solid, pools words
weighting by how much each matters to the baseline — frequent words dominate)
and **w-Macro** (dotted, every word counts equally), with their gap shaded.

**How to read it.** A large solid-over-dotted gap means the model fits **common
words much better than rare ones**. This is the exact word-space analogue of
E3's length-bias gap: there, long documents were fit better; here, frequent
words are.

**Takeaway.** w-Micro sits far above w-Macro at every K (at the full design
point: 0.34 vs 0.19 at K\* = 40, gap ≈ 0.15), and both curves elbow at K\*. A
persistent w-Micro ≫ w-Macro gap is the signal that "good overall fit" is
really "good fit on frequent words."

### `F12_word_fit_vs_freq` — where in the frequency spectrum fit lives
**What you see.** Per-word held-out R² at K\* against log document-frequency,
with a trend line.

**How to read it.** The upward trend confirms frequent words fit better;
points far below the trend are the individually mispredicted words worth
inspecting (the word-level analogue of E4's `F5b` residual ranking).

**Takeaway.** Fit concentrates in the common vocabulary; the rare tail is
estimation-hard. A companion table also verifies **Lemma 2** — the total
document-wise and word-wise deviance agree to ~10⁻⁷ or better, confirming the
two views are the same fit seen through different weights.

---

## Figures that appear once their runs complete

- **`F10` (E1b/E1c)** — the E1 story repeated at higher complexity (`K* = 20`) and with sparser documents (`α_DGP = 0.1`): does the elbow stay sharp? The variants are **manual-only** now (`Rscript Code/run_E1.R full E1b <workers>`); the figure populates when they run.

(`F8` no longer exists — see its entry above.)

---

*Generated for the pilot run. Regenerate figures with
`Rscript Code/make_figures.R <profile>`; see `Code/README.md` for the full
figure/experiment catalog and `Results/PILOT_REPORT.md` for the numeric summary.*
