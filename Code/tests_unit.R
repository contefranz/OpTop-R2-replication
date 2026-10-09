# =============================================================================
# tests_unit.R -- correctness gates for the Section 5 harness.
# Run BEFORE any experiment:  Rscript Code/tests_unit.R
# Exits non-zero on first failure.
#
#   U1 local held-out partition == OpTop::optop_make_partition (in-sample)
#   U2 pseudo-fit scoring path == native OpTop path (indices, 1e-10)
#   U3 Prop-1(ii) identity: influence-function gap == Micro - Macro
#   U4 Prop-1(iii) decomposition sums to the gap at machine precision
#   U5 token split conserves counts; fold-in side non-empty for L >= 2
#   U6 Definition-1 selector picks the documented K on a toy gains table
#      (and is not derailed by an NA bound elsewhere on the grid)
#   U7 instruments: rows sum to zero; moment test ~ nominal size on iid noise
#   U8 EM fold-in ~ topicmodels posterior        U9 WarpLDA slim-cache path
#   U10 Lemma S1: PACKAGE word-level deviances == direct document-wise total
#   U11 word over-K helper + planted words       U12-U13 DGP vs fit priors
#   U14 fit-cache key invariance                 U15 native == reference partition
#   U16 null-discrepancy floor
#   U17 Test-3 instruments are K-specific (data.table shadowing regression)
#   U18 simultaneous / total-gain selectors; cluster SE; "none certified"
#   U19 held-out word-level null deviance in Poisson form (Section 2.4)
# =============================================================================

source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

fail <- function(...) { cat("FAIL:", sprintf(...), "\n"); quit(status = 1L) }
pass <- function(id) cat(sprintf("PASS %s\n", id))

set.seed(1)
sim <- sim_lda_corpus(J = 80L, W = 400L, K_true = 4L, alpha = 0.5, beta = 0.05,
                      length_spec = list(type = "fixed", L = 200L), seed = 42L)
future::plan(future::sequential)
fits <- fit_lda_grid(sim$dtm, 2:5, method = "VEM", n_starts = 1L,
                     fit_seed_base = 7L, parallel = FALSE)$models

# --- U1/U2: cross-check against the package ------------------------------------
chk <- crosscheck_optop(fits, sim$dtm, c = 1, tol = 1e-10)
if (!all(unlist(chk))) fail("crosscheck_optop: %s",
                            paste(names(chk)[!unlist(chk)], collapse = ", "))
pass("U1/U2 partition + pseudo-fit path == native OpTop")

# --- U3/U4: gap identity and decomposition ---------------------------------------
ins <- score_insample(fits, sim$dtm, c = 1, metrics = "dev")
gap_if <- micro_macro_gap_ci(ins$doc)
cmp <- merge(gap_if, ins$summary[, .(K, metric, gap_sum = r2_micro - r2_macro)],
             by = c("K", "metric"))
if (cmp[, max(abs(gap - gap_sum))] > 1e-10)
  fail("U3 gap identity: max dev %.3e", cmp[, max(abs(gap - gap_sum))])
pass("U3 influence-function gap == Micro - Macro")

dec <- gap_decomposition(ins$doc[metric == "dev"])
if (dec[, max(abs(resid))] > 1e-12)
  fail("U4 decomposition residual %.3e", dec[, max(abs(resid))])
pass("U4 Prop-1(iii) channels sum to gap")

# --- U5: token split ---------------------------------------------------------------
spl <- split_tokens_binomial(sim$dtm, 0.5, seed = 3L)
if (!identical(as.matrix(spl$foldin + spl$score), as.matrix(sim$dtm)))
  fail("U5 split does not conserve counts")
L <- Matrix::rowSums(sim$dtm)
if (any(Matrix::rowSums(spl$foldin)[L >= 2] == 0))
  fail("U5 empty fold-in side for a document with L >= 2")
pass("U5 token split conserves counts, fold-in non-empty")

# --- U6: Definition-1 selector -------------------------------------------------------
toy <- data.table(metric = "dev", K = 2:6, K_next = 3:7,
                  delta_mean = c(.20, .10, .004, .002, .001),
                  delta_sd = .01, n = 100)
toy[, se := delta_sd / sqrt(n)]
toy[, `:=`(ub_onesided = delta_mean + qnorm(.95) * se, alpha = 0.05)]
kh <- select_k_epsilon(toy, eps = 0.01, alpha = 0.05)
if (kh$K_hat != 4L) fail("U6 expected K_hat = 4, got %s", kh$K_hat)
# an undefined bound at a LATER grid point must not erase a valid selection
toy_na <- copy(toy)[K == 6L, ub_onesided := NA_real_]
kh_na <- select_k_epsilon(toy_na, eps = 0.01, alpha = 0.05)
if (!identical(kh_na$K_hat, 4L))
  fail("U6 NA bound at K = 6 changed the selection to %s", kh_na$K_hat)
pass("U6 epsilon-rule selects the documented K (robust to NA bounds)")

# --- U7: instruments and Wald size on iid noise ---------------------------------------
Z1 <- make_instrument_freq_contrast(sim$dtm)
Z2 <- make_instruments_freq_strata(sim$dtm, B = 5L)
if (max(abs(rowSums(Z1)), abs(rowSums(Z2))) > 1e-12)
  fail("U7 instrument rows do not sum to zero")
set.seed(9)
pvals <- replicate(300, {
  E <- matrix(rnorm(150 * ncol(sim$dtm), sd = 1e-3), nrow = 150)
  E <- E - rowMeans(E)                       # residuals sum to zero over vocab
  moment_test(Z2, E)$result$pval
})
sz <- mean(pvals < 0.05)
if (sz < 0.01 || sz > 0.12) fail("U7 null size %.3f outside [0.01, 0.12]", sz)
pass(sprintf("U7 instruments zero-sum; Wald null size %.3f", sz))

# --- U8: EM fold-in agrees with topicmodels posterior ----------------------------
fitv <- fits[["4"]]
th_post <- foldin_theta(fitv, sim$dtm)
th_em <- foldin_theta_em(phi_from_fit(fitv), sim$dtm)
cc <- cor(as.vector(th_post), as.vector(th_em))
if (cc < 0.98) fail("U8 EM vs posterior fold-in corr %.4f < 0.98", cc)
pass(sprintf("U8 EM fold-in ~ posterior (corr %.4f; ML vs Dirichlet-prior pull)", cc))

# --- U9: WarpLDA end-to-end through slim cache -------------------------------------
fw <- fit_lda_grid(sim$dtm, 2:5, method = "WarpLDA", n_starts = 2L,
                   fit_seed_base = 7L, parallel = FALSE)$models
if (!all(vapply(fw, function(f) is.null(f$model_object), TRUE)))
  fail("U9 WarpLDA fits were not slimmed before caching")
tf <- tempfile(fileext = ".qs2")
qs2::qs_save(fw, tf); fw2 <- qs2::qs_read(tf)          # cache round-trip
insw <- score_insample(fw2, sim$dtm, c = 1, metrics = "dev")
if (!all(is.finite(insw$summary$r2_micro)))
  fail("U9 WarpLDA in-sample indices not finite")
basew <- OpTop::optop_make_baseline(sim$dtm)
how <- score_heldout(fw2, sim$dtm, sim$dtm, basew$pi_glob, c = 1, "dev")
decw <- gap_decomposition(how$doc[metric == "dev"])
if (decw[, max(abs(resid))] > 1e-12)
  fail("U9 WarpLDA gap decomposition residual %.3e", decw[, max(abs(resid))])
gw <- paired_gains(how$doc, alpha = 0.05)
khw <- select_k_epsilon(gw, eps = 0.05, alpha = 0.05)
if (!is.integer(khw$K_hat)) fail("U9 K-hat machinery failed on WarpLDA path")
pass("U9 WarpLDA: slim-cache round-trip, EM fold-in, indices + K-hat machinery")

# --- U10: Lemma S1 (doc-wise == word-wise unbinned fitted deviance) ----------------
# The identity is checked ACROSS CODE PATHS: a directly computed document-wise
# total against the package's own word-level deviances, in total and word by
# word. (Re-adding one vector of cell terms in two orders -- the pre-revision
# gate -- agrees by construction and could never fail.)
sc10 <- score_insample(fits, sim$dtm, c = 1, metrics = "dev",
                       word_at = 4L, word_metrics = "dev")
w10 <- sc10$word[K == 4L & metric == "dev"]
w10 <- w10[match(colnames(sim$dtm), word_id)]
if (anyNA(w10$word_id)) fail("U10 word table does not cover the vocabulary")
l2 <- lemma2_residual(theta_from_fit(fits[["4"]]), phi_from_fit(fits[["4"]]),
                      sim$dtm, word_d_model = w10$d_model)
if (abs(l2$rel_resid_pkg) > 1e-8)
  fail("U10 doc-wise total %.6f vs package word-wise total %.6f (rel %.3e)",
       l2$dev_doc, l2$dev_word_pkg, l2$rel_resid_pkg)
if (l2$max_word_diff > 1e-6 * max(1, max(abs(w10$d_model))))
  fail("U10 package word deviance differs from the direct Poisson form (max %.3e)",
       l2$max_word_diff)
pass(sprintf("U10 Lemma S1 across code paths (rel. total %.1e; max per-word %.1e)",
             l2$rel_resid_pkg, l2$max_word_diff))

# --- U11: word-level over-K helper + planted-word detection ------------------------
# High-frequency stopword contamination (E6 scenario W2): the planted words are
# non-topic common words the model predicts no better than the global baseline,
# so they populate the worst-fit tail.
Phi_h <- sim$Phi
pi_stop_h <- make_stopword_dist(400L, n_stop = 20L, seed = 5L)
planted_h <- which(pi_stop_h > 0)
sim_h <- sim_lda_corpus(120L, 400L, 4L, 0.5, 0.05,
                        list(type = "fixed", L = 250L), seed = 43L,
                        Phi = Phi_h, doc_prefix = "d",
                        contamination = list(pi_stop = pi_stop_h, w_mean = 0.2,
                                             w_conc = 10, mode = "shared"))
fits_h <- fit_lda_grid(sim_h$dtm, 3:5, method = "VEM", n_starts = 1L,
                       fit_seed_base = 11L, parallel = FALSE)$models
base_h <- OpTop::optop_make_baseline(sim_h$dtm)
sc_h <- score_insample(fits_h, sim_h$dtm, c = 1, metrics = "dev",
                       word_at = 3:5, word_metrics = "dev")
wc <- word_micro_macro_over_k(sc_h$word, sim_h$dtm, base_h$pi_glob, 3L, 1)$curve
if (!nrow(wc) || any(!is.finite(wc$gap)) || any(wc$n_words <= 0))
  fail("U11 word_micro_macro_over_k returned no finite indices")
wstar <- word_index_filtered(sc_h$word[K == 4L], sim_h$dtm, base_h$pi_glob,
                             3L, 1)$word[keep == TRUE][order(r2_word)]
is_planted <- wstar$word_id %in% paste0("word_", planted_h)
worst_third <- mean(is_planted[seq_len(ceiling(nrow(wstar) / 3))])
if (!(worst_third > mean(is_planted)))
  fail("U11 planted words not enriched in the worst-fit tail (%.2f vs %.2f)",
       worst_third, mean(is_planted))
pass(sprintf("U11 word over-K helper finite; planted enriched in worst tail (%.2f vs %.2f)",
             worst_third, mean(is_planted)))

# --- U12: DGP priors and fitted-model priors stay distinct -----------------------
# The simulation config deliberately uses alpha_DGP/beta_DGP, while alpha/beta
# are optional fitting priors (NA means engine default).  This catches a silent
# but severe failure mode: accidentally feeding NA fitting priors into the DGP.
cfgp <- get_config("E1", "smoke")
sim_p <- sim_lda_corpus(20L, 160L, cfgp$K_true,
                        cfgp$alpha_DGP, cfgp$beta_DGP,
                        list(type = "fixed", L = 100L), seed = 99L)
if (!all(is.finite(sim_p$Phi)) ||
    max(abs(rowSums(sim_p$Phi) - 1)) > 1e-12 || sum(sim_p$Phi > 0) <= cfgp$K_true)
  fail("U12 DGP priors did not produce a non-degenerate topic matrix")

fit_p <- fit_lda_grid(sim$dtm, 4L, method = "VEM", n_starts = 1L,
                      fit_seed_base = 17L, fit_alpha = 0.3, parallel = FALSE)
alpha_p <- mean(fit_p$models[["4"]]$model_object@alpha)
if (abs(alpha_p - 0.3) > 1e-10)
  fail("U12 fixed VEM alpha %.4g, expected 0.3", alpha_p)
bad_beta <- try(fit_lda_grid(sim$dtm, 4L, method = "VEM", n_starts = 1L,
                             fit_seed_base = 17L, fit_beta = 0.01,
                             parallel = FALSE), silent = TRUE)
if (!inherits(bad_beta, "try-error"))
  fail("U12 VEM accepted an unsupported fixed beta")

tag_base <- run_tag(cfgp)
cfgp_dgp <- cfgp; cfgp_dgp$alpha_DGP <- 0.1
cfgp_fit <- cfgp; cfgp_fit$alpha <- 0.1
if (identical(tag_base, run_tag(cfgp_dgp)) || identical(tag_base, run_tag(cfgp_fit)) ||
    identical(run_tag(cfgp_dgp), run_tag(cfgp_fit)))
  fail("U12 run tags do not distinguish DGP and fitted-model priors")
pass("U12 DGP/fit priors are distinct, VEM alpha is fixed, VEM beta fails fast")

# --- U13: engines that expose both fitted priors receive both --------------------
fit_g <- fit_lda_grid(sim$dtm, 4L, method = "Gibbs", n_starts = 1L,
                      fit_seed_base = 23L,
                      gibbs = list(burnin = 10L, iter = 30L, thin = 5L),
                      fit_alpha = 0.4, fit_beta = 0.02, parallel = FALSE)
if (abs(mean(fit_g$models[["4"]]$model_object@alpha) - 0.4) > 1e-10 ||
    abs(fit_g$models[["4"]]$model_object@control@delta - 0.02) > 1e-10)
  fail("U13 Gibbs did not receive the configured alpha/delta priors")
fit_w <- fit_lda_grid(sim$dtm, 4L, method = "WarpLDA", n_starts = 1L,
                      fit_seed_base = 23L,
                      fit_alpha = 0.4, fit_beta = 0.02, parallel = FALSE)
if (fit_w$diagnostics$fit_alpha != 0.4 || fit_w$diagnostics$fit_beta != 0.02)
  fail("U13 WarpLDA did not retain the configured priors")
# Behavioral check: the diagnostics columns above merely echo the inputs, so
# also prove the prior REACHES the sampler -- same seed, very different
# doc_topic_prior, materially different fitted theta.
fit_w2 <- fit_lda_grid(sim$dtm, 4L, method = "WarpLDA", n_starts = 1L,
                       fit_seed_base = 23L,
                       fit_alpha = 5.0, fit_beta = 0.02, parallel = FALSE)
th_lo <- theta_from_fit(fit_w$models[["4"]])
th_hi <- theta_from_fit(fit_w2$models[["4"]])
if (max(abs(th_lo - th_hi)) < 1e-6)
  fail("U13 WarpLDA theta insensitive to doc_topic_prior (0.4 vs 5.0) -- prior not reaching the engine")
pass("U13 Gibbs and WarpLDA receive explicit fitted priors (WarpLDA verified behaviorally)")

# --- U14: fit-cache key invariance ------------------------------------------------
# (a) With fit priors unset, dgp_signature() must keep the PRE-RENAME field
#     layout byte-for-byte, or every legacy Data/FITS hash silently rots.
# (b) Two independently built identical configs -> the same cache path.
# (c) Setting a fit prior -> a different cache path (no stale-fit reuse).
seeds14 <- make_seeds(cfgp$seed_base, 1L)
sig_unset <- dgp_signature(cfgp, seeds14)
if (!identical(names(sig_unset),
               c("K_true", "J", "W", "alpha", "beta", "length_spec",
                 "dgp_seed", "scenario", "extra")))
  fail("U14 dgp_signature layout changed for unset priors (legacy FITS hashes would rot): %s",
       paste(names(sig_unset), collapse = ", "))
p14_a <- fit_cache_path_k(sig_unset, 4L, "VEM", 1L, 17L)
p14_b <- fit_cache_path_k(dgp_signature(get_config("E1", "smoke"), seeds14),
                          4L, "VEM", 1L, 17L)
if (!identical(p14_a, p14_b))
  fail("U14 identical configs give different fit-cache paths")
cfgp14 <- cfgp; cfgp14$alpha <- 0.3
p14_c <- fit_cache_path_k(dgp_signature(cfgp14, seeds14), 4L, "VEM", 1L, 17L)
if (identical(p14_a, p14_c))
  fail("U14 setting a fit prior did not change the fit-cache path")
pass("U14 fit-cache key: legacy layout preserved when unset; prior changes the key")

# --- U15: OpTop 0.14.0 native partition == pure-R reference; holdout agreement ---
# The production make_heldout_partition() delegates to
# OpTop::optop_make_partition(pi_glob=) since 0.14.0; the pre-0.14.0 pure-R
# builder is retained as .make_heldout_partition_ref(). (a) They must agree
# exactly on a HELD-OUT case (eval dtm scored under the TRAINING baseline).
# (b) The full local scorer must agree with the package's own
# optop_index_holdout() on live fits (document-level reconstruction, no OOV).
sim_ev15 <- sim_lda_corpus(J = 40L, W = 400L, K_true = 4L, alpha = 0.5,
                           beta = 0.05,
                           length_spec = list(type = "fixed", L = 200L),
                           seed = 43L, Phi = sim$Phi)
base15 <- OpTop::optop_make_baseline(sim$dtm)
th15 <- lapply(fits, foldin_theta, newdata_dtm = sim_ev15$dtm, seed = 5L)
ph15 <- lapply(fits, phi_from_fit)
p_nat <- make_heldout_partition(th15, ph15, sim_ev15$dtm, base15$pi_glob, c = 1)
p_ref <- .make_heldout_partition_ref(th15, ph15, sim_ev15$dtm, base15$pi_glob,
                                     c = 1)
if (!identical(unname(p_nat$rare_mask) > 0, unname(p_ref$rare_mask) > 0) ||
    !isTRUE(all.equal(as.numeric(p_nat$L), as.numeric(p_ref$L),
                      tolerance = 1e-10)) ||
    !identical(as.logical(p_nat$chisq_min_ok), as.logical(p_ref$chisq_min_ok)) ||
    !isTRUE(all.equal(p_nat$chisq_min_report$excluded_mass,
                      p_ref$chisq_min_report$excluded_mass, tolerance = 1e-10)))
  fail("U15 native held-out partition differs from the pure-R reference")
ho15 <- OpTop::optop_index_holdout(models = unname(fits),
                                   dtm_eval = sim_ev15$dtm,
                                   baseline = base15, c = 1,
                                   metrics = "deviance")
loc15 <- score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                       c = 1, metrics = "dev")
nat15 <- as.data.frame(ho15$summary)[order(ho15$summary$K), ]
lcl15 <- as.data.frame(loc15$summary)[order(loc15$summary$K), ]
if (!isTRUE(all.equal(nat15$micro, lcl15$r2_micro, tolerance = 1e-8)) ||
    !isTRUE(all.equal(nat15$macro, lcl15$r2_macro, tolerance = 1e-8)))
  fail("U15 local held-out scorer disagrees with optop_index_holdout")
pass("U15 native partition == reference; local scorer == optop_index_holdout")

# --- U16: OpTop 0.14.1 null-discrepancy floor -------------------------------------
# (a) Direct engine semantics on a 3-doc synthetic: healthy / soft-degenerate /
#     exact-collapse. min_null = 1 excludes the two degenerate docs from BOTH
#     aggregations; min_null = 0 reproduces the legacy strict-positivity rule.
D_K16    <- c(5,    0.4,  1e-21)
D_null16 <- c(20,   1e-5, 4e-27)
r16 <- OpTop:::.optop_index_result_doc(D_K16, D_null16, 4L, "deviance",
                                       macro = TRUE,   # ztest dropped in 0.20
                                       min_null = 1)
if (r16$n_null_excluded != 2L ||
    abs(r16$r2_macro - (1 - 5 / 20)) > 1e-12 ||
    !identical(is.na(r16$r2_doc), c(FALSE, TRUE, TRUE)))
  fail("U16 null floor did not exclude degenerate docs / wrong Macro")
r16L <- OpTop:::.optop_index_result_doc(D_K16, D_null16, 4L, "deviance",
                                        macro = TRUE,  # ztest dropped in 0.20
                                        min_null = 0)
if (r16L$n_null_excluded != 0L || any(is.na(r16L$r2_doc)))
  fail("U16 min_null = 0 did not reproduce the legacy behavior")
# (b) Pipeline path: healthy toy corpus -> the floor is a no-op and the summary
#     reports null_excl_share = 0; blocked and unblocked pooling agree.
sc16 <- score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                      c = 1, metrics = "dev")
if (!"null_excl_share" %in% names(sc16$summary) ||
    any(sc16$summary$null_excl_share > 0))
  fail("U16 healthy corpus reported a nonzero null-floor exclusion share")
bl16 <- score_heldout_blocked(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                              c = 1, metrics = "dev", block_docs = 20L)
cmp16 <- merge(sc16$summary[, .(K, r2_micro, r2_macro)],
               bl16$summary[, .(K, r2_micro_b = r2_micro, r2_macro_b = r2_macro)],
               by = "K")
if (!isTRUE(all.equal(cmp16$r2_micro, cmp16$r2_micro_b, tolerance = 1e-10)) ||
    !isTRUE(all.equal(cmp16$r2_macro, cmp16$r2_macro_b, tolerance = 1e-10)))
  fail("U16 blocked pooling disagrees with the floored engine")
pass("U16 null-discrepancy floor: excludes degenerate docs, legacy at 0, blocked == unblocked")

# --- U17: Test-3 instruments are K-specific ----------------------------------------
# Regression gate for a data.table scoping defect: inside `[.data.table` a
# closure argument named like a column is shadowed by the column, so
# `word[K == K & ...]` keeps every K and `lookup[list(K), ws][[1L]]` self-joins;
# both silently gave every tested K the strata of the SMALLEST one. The gate
# exercises the shared builder the drivers call, not a copy of it.
sc17 <- score_insample(fits, sim$dtm, c = 1, metrics = "dev",
                       word_at = 2:5, word_metrics = "dev")
Kt17 <- c(2L, 4L, 5L)
Z17 <- make_instruments_by_K(sim$dtm, sc17$word, Kt17, B = 5L, S = 5L,
                             min_docfreq = 2L)
df17 <- as.integer(Matrix::colSums(sim$dtm > 0))
for (k17 in Kt17) {
  ref17 <- make_instruments_fit_strata(
    sc17$word[K == k17 & metric == "dev", .(word_id, r2_word)],
    colnames(sim$dtm), S = 5L, min_docfreq = 2L, docfreq = df17)
  if (!identical(Z17[[as.character(k17)]]$T3_fit_strata, ref17))
    fail("U17 Test-3 instruments at K = %d are not built from that K's scores", k17)
  if (max(abs(rowSums(Z17[[as.character(k17)]]$T3_fit_strata))) > 1e-12)
    fail("U17 Test-3 instrument rows do not sum to zero at K = %d", k17)
}
if (identical(Z17[["2"]]$T3_fit_strata, Z17[["5"]]$T3_fit_strata))
  fail("U17 Test-3 instruments identical at K = 2 and K = 5 (shadowing?)")
if (!identical(Z17[["2"]]$T2_freq_strata, Z17[["5"]]$T2_freq_strata))
  fail("U17 frequency instruments should not depend on K")
bad17 <- try(make_instruments_by_K(sim$dtm, sc17$word, 99L), silent = TRUE)
if (!inherits(bad17, "try-error"))
  fail("U17 a K without training word scores did not fail loudly")
pass("U17 Test-3 instruments are K-specific; frequency instruments K-invariant")

# --- U18: simultaneous and total-gain selectors --------------------------------------
# Plateau with a late step: adjacent gains .2, .1, .004, .001, so the gain still
# available from K = 4 is .005. The LOCAL rule certifies 4 at eps = .0045 (next
# step .004) while TOTAL-GAIN cannot (.005 remains) and moves to 5 -- the
# distinction the revised Definition 1 draws.
set.seed(18)
n18 <- 200L; K18 <- 2:6; base18 <- c(.30, .50, .60, .604, .605)
u18 <- rnorm(n18, sd = .05)
doc18 <- rbindlist(lapply(seq_along(K18), function(i) data.table(
  K = K18[i], metric = "dev", doc_id = sprintf("d%03d", seq_len(n18)),
  r2_doc = base18[i] + u18 + rnorm(n18, sd = 1e-4))))
pa18 <- paired_gains_all(doc18, alpha = 0.05, family = "all_pairs")
ad18 <- paired_gains_all(doc18, alpha = 0.05, family = "adjacent")
if (pa18$M[1L] != 10L || ad18$M[1L] != 4L || nrow(pa18) != 10L)
  fail("U18 comparison-family sizes wrong (M = %d / %d)", pa18$M[1L], ad18$M[1L])
tg_a <- select_k_total_gain(pa18, eps = 0.01)
tg_b <- select_k_total_gain(pa18, eps = 0.0045)
tg_c <- select_k_total_gain(pa18, eps = 1e-5)
lo_b <- select_k_adjacent(ad18, eps = 0.0045, simultaneous = TRUE)
if (!identical(tg_a$K_hat, 4L) || tg_a$n_comparisons != 2L)
  fail("U18 total-gain at eps = .01: K_hat %s, comparisons %s",
       tg_a$K_hat, tg_a$n_comparisons)
if (!identical(tg_b$K_hat, 5L) || tg_b$n_comparisons != 1L ||
    !identical(lo_b$K_hat, 4L))
  fail("U18 local/total distinction lost (total %s, local %s)",
       tg_b$K_hat, lo_b$K_hat)
if (!is.na(tg_c$K_hat) || isTRUE(tg_c$certified))
  fail("U18 empty certified set must give NA / certified = FALSE, never the boundary")
if (any(pa18$K == max(K18)))
  fail("U18 the grid maximum must never be a candidate")
# cached-table path == document-table path
lo_cached <- select_k_adjacent_simultaneous(paired_gains(doc18, alpha = 0.05),
                                            eps = 0.0045, alpha = 0.05)
if (!identical(lo_cached$K_hat, lo_b$K_hat) ||
    abs(lo_cached$max_ub - lo_b$max_ub) > 1e-12)
  fail("U18 cached-gains selector disagrees with the document-table selector")
# singleton clusters: cluster-robust SE == iid SE exactly
cl18 <- data.table(doc_id = unique(doc18$doc_id))[, cluster_id := doc_id]
pc18 <- paired_gains_all(doc18, alpha = 0.05, cluster = cl18, family = "all_pairs")
if (max(abs(pc18$se - pa18$se)) > 1e-14 || pc18$se_type[1L] != "cluster")
  fail("U18 singleton-cluster SE differs from the iid SE (max %.3e)",
       max(abs(pc18$se - pa18$se)))
# perfectly duplicated documents within cluster: SE must grow by sqrt(2)-ish
cl18b <- copy(cl18)[, cluster_id := sprintf("g%03d", (seq_len(.N) + 1L) %/% 2L)]
doc18b <- copy(doc18)[, pair := (as.integer(sub("d", "", doc_id)) + 1L) %/% 2L]
doc18b[, r2_doc := mean(r2_doc), by = .(K, pair)][, pair := NULL]
se_i <- paired_gains_all(doc18b, family = "adjacent")$se
se_c <- paired_gains_all(doc18b, cluster = cl18b, family = "adjacent")$se
if (any(se_c / se_i < 1.35 | se_c / se_i > 1.48))
  fail("U18 cluster SE ratio under duplicated pairs %.3f, expected ~sqrt(2)",
       mean(se_c / se_i))
# a retained set that varies with K is rejected
bad18 <- copy(doc18)[K == 4L & doc_id == "d001", r2_doc := NA_real_]
if (!inherits(try(paired_gains_all(bad18), silent = TRUE), "try-error"))
  fail("U18 K-dependent retained set was not rejected")
ar18 <- select_k_all_rules(doc18, eps_grid = c(0.01, 0.0045))
if (nrow(ar18) != 6L || !all(c("adjacent_pointwise", "adjacent_simultaneous",
                                "total_gain") %in% ar18$rule))
  fail("U18 select_k_all_rules layout")
pass("U18 selectors: family sizes, local vs total gain, none-certified, cluster SE")

# --- U19: held-out word-level baseline deviance is in Poisson form ------------------
# Section 2.4 defines D_w(null) = 2 sum_j [N log(N/B) - (N - B)]. OpTop <= 0.20.1
# omits the linear term on the word-level NULL path; in-sample that term is
# identically zero, held-out (training baseline, evaluation counts) it is not.
# The pipeline therefore computes the held-out null itself; the gate checks it
# against a dense, direct evaluation of the definition.
N19 <- as.matrix(sim_ev15$dtm); L19 <- rowSums(N19)
B19 <- outer(L19, as.numeric(base15$pi_glob))
ok19 <- colSums(N19) > 0 & as.numeric(base15$pi_glob) > 0   # no flooring involved
T19 <- ifelse(N19 > 0, N19 * log(N19 / B19), 0) - (N19 - B19)
ref19 <- 2 * colSums(T19)
loc19 <- .word_null_dev_poisson(sim_ev15$dtm, base15$pi_glob)
if (max(abs(loc19[ok19] - ref19[ok19])) > 1e-8 * max(1, max(abs(ref19[ok19]))))
  fail("U19 held-out word null deviance differs from the Poisson definition (max %.3e)",
       max(abs(loc19[ok19] - ref19[ok19])))
if (min(ref19[ok19]) < -1e-8) fail("U19 Poisson null deviance must be nonnegative")
w19 <- score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob, c = 1,
                     metrics = "dev", word_at = 4L)$word
w19 <- w19[match(colnames(sim_ev15$dtm), word_id)]
if (max(abs(w19$d_null[ok19] - ref19[ok19])) > 1e-8 * max(1, max(abs(ref19[ok19]))))
  fail("U19 score_heldout() word table does not carry the Poisson-form null")
chk19 <- ok19 & w19$d_null > 0
if (max(abs(w19$r2_word[chk19] - (1 - w19$d_model[chk19] / w19$d_null[chk19]))) > 1e-12)
  fail("U19 r2_word is not 1 - d_model / d_null")
# in-sample the linear term vanishes: the package null is already the Poisson form
ins19 <- score_insample(fits, sim$dtm, c = 1, metrics = "dev", word_at = 4L)$word
ins19 <- ins19[match(colnames(sim$dtm), word_id)]
loc19i <- .word_null_dev_poisson(sim$dtm, OpTop::optop_make_baseline(sim$dtm)$pi_glob)
pos19 <- Matrix::colSums(sim$dtm) > 0
if (max(abs(ins19$d_null[pos19] - loc19i[pos19])) > 1e-8 * max(1, max(abs(loc19i))))
  fail("U19 in-sample: package null and Poisson-form null should coincide")
pass("U19 held-out word-level null deviance in Poisson form; coincides in-sample")


# --- U20: the package's validity warning, and the correction applied exactly once ----
# OpTop 0.20.1 signals a classed warning (optop_word_null_baseline) when a
# word-level null is computed on an EXTERNAL baseline; its kernel is unchanged.
# (a) the raw package call warns; (b) score_heldout() replaces that null by the
# Poisson form and muffles THIS warning only, in the branch that replaces it;
# (c) when the fitted probabilities equal the baseline every index is zero, under
# reconstruction and completion alike; (d) converting legacy -> Poisson and then
# passing the result through again changes nothing (no double correction).
N20 <- as(as(Matrix(matrix(c(3, 0, 0, 2), 2, 2, byrow = TRUE,
  dimnames = list(c("d1", "d2"), c("a", "b"))), sparse = TRUE),
  "generalMatrix"), "CsparseMatrix")
p20 <- c(a = .2, b = .8)
f20 <- as_pseudo_fit(matrix(1, 2, 1), matrix(p20, 1, 2), rownames(N20), colnames(N20))
for (target20 in c("rec", "com")) {
  S20 <- N20
  if (target20 == "com") S20@x[] <- 1          # scored tokens differ from the document
  part20 <- make_heldout_partition(list(matrix(1, 2, 1)), list(matrix(p20, 1, 2)),
                                   S20, p20, .1)
  caught20 <- FALSE
  withCallingHandlers(
    OpTop::optop_index_deviance(f20, S20, part20, list(pi_glob = p20), level = "word"),
    optop_word_null_baseline = function(w) {
      caught20 <<- TRUE; invokeRestart("muffleWarning") })
  if (!caught20) fail("U20 OpTop did not signal optop_word_null_baseline (%s)", target20)
  out20 <- withCallingHandlers(
    score_heldout(list("1" = f20), N20, S20, p20, c = .1, metrics = "dev",
                  word_at = 1L, min_null = 0),
    optop_word_null_baseline = function(w)
      fail("U20 the package warning leaked out of the corrected branch"))
  if (!identical(out20$scoring$word_null_convention, "poisson_scored_tokens"))
    fail("U20 score_heldout() result is not stamped with its word-null convention")
  B20 <- outer(as.numeric(rowSums(S20)), p20); A20 <- as.matrix(S20)
  ref20 <- 2 * colSums(ifelse(A20 > 0, A20 * log(pmax(A20, 1) / B20), 0) - A20 + B20)
  if (max(abs(out20$word$d_null - ref20)) > 1e-10 ||
      max(abs(out20$word$r2_word)) > 1e-10 || max(abs(out20$doc$r2_doc)) > 1e-10)
    fail("U20 fitted == baseline must give a zero index on the Poisson null (%s)", target20)
}
# an UNRELATED warning raised inside the same call must survive the handler
leak20 <- FALSE
withCallingHandlers(
  withCallingHandlers({ warning("unrelated"); 1 },
    optop_word_null_baseline = function(w) invokeRestart("muffleWarning")),
  warning = function(w) { leak20 <<- TRUE; invokeRestart("muffleWarning") })
if (!leak20) fail("U20 a classed handler must not swallow other warnings")
Nw20 <- as.numeric(colSums(S20)); Bw20 <- p20 * sum(S20)
legacy20 <- copy(out20$word)[, d_null := d_null + 2 * (Nw20 - Bw20)]
a20 <- revision_word_null(legacy20, Nw20, Bw20, "legacy_log_only")
b20 <- revision_word_null(a20, Nw20, Bw20, "poisson_scored_tokens")
if (!identical(a20, b20) || max(abs(a20$d_null - ref20)) > 1e-10)
  fail("U20 the legacy correction is not idempotent through the convention guard")
if (!inherits(try(revision_word_null(a20, Nw20, Bw20, "unknown"), silent = TRUE), "try-error") ||
    !inherits(try(revision_word_null(a20, Nw20, Bw20), silent = TRUE), "try-error"))
  fail("U20 an unknown or missing convention must be an error")
pass("U20 package guard caught and confined; zero index at fitted == baseline (rec, com); correction exactly once")

# --- U21: c (support) and delta (floor) are independent ---------------------------------
lo21 <- score_heldout(list("1" = f20), N20, N20, p20, c = .1, metrics = "dev", min_null = 0)
hi21 <- suppressMessages(
  score_heldout(list("1" = f20), N20, N20, p20, c = .1, metrics = "dev", min_null = 2))
if (!identical(lo21$partition, hi21$partition) ||
    !identical(lo21$doc$d_null, hi21$doc$d_null) ||
    !identical(lo21$doc$d_model, hi21$doc$d_model))
  fail("U21 the discrepancy floor changed the support or a raw discrepancy")
if (!(sum(is.na(hi21$doc$r2_doc)) > sum(is.na(lo21$doc$r2_doc))))
  fail("U21 a higher floor must exclude more documents")
# ... and on the shared fixture: moving c moves the support, the floor does not
sc21a <- suppressMessages(score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                                        c = 0.5, metrics = "dev", min_null = 1))
sc21b <- suppressMessages(score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                                        c = 2, metrics = "dev", min_null = 1))
sc21c <- suppressMessages(score_heldout(fits, sim_ev15$dtm, sim_ev15$dtm, base15$pi_glob,
                                        c = 2, metrics = "dev", min_null = 3))
if (identical(sc21a$partition$rare_mask, sc21b$partition$rare_mask))
  fail("U21 changing c did not change the support")
if (!identical(sc21b$partition$rare_mask, sc21c$partition$rare_mask) ||
    !identical(sc21b$doc$d_null, sc21c$doc$d_null))
  fail("U21 changing delta at fixed c changed the support or the discrepancies")
pass("U21 support threshold c and discrepancy floor delta act independently")

# --- U22: checkpoints are reused only on identical configuration / package / code ------
old22 <- getOption("optop.revision_suffix", "")
sfx22 <- paste0("_unit", Sys.getpid())
options(optop.revision_suffix = sfx22)
v22 <- revision_checkpoint("unit", list(seed = 22), "one", 42L)
if (!identical(revision_checkpoint("unit", list(seed = 22), "one", stop("recomputed")), 42L))
  fail("U22 a stored checkpoint was recomputed")
ck22 <- list.files(p_data("Checkpoints", sfx22), "\\.qs2$", recursive = TRUE, full.names = TRUE)
tam22 <- qs2::qs_read(ck22[1L]); tam22$identity$code <- "another version of the scoring code"
qs2::qs_save(tam22, ck22[1L])
if (!inherits(try(revision_checkpoint("unit", list(seed = 22), "one", 42L), silent = TRUE),
              "try-error"))
  fail("U22 a checkpoint written by different code must be refused, not reused")
unlink(p_data("Checkpoints", sfx22), recursive = TRUE)
options(optop.revision_suffix = "")
if (!identical(revision_checkpoint("unit", list(seed = 22), "one", 7L), 7L) ||
    dir.exists(file.path(proj_path("Data", "Checkpoints"), sfx22)))
  fail("U22 without a suffix nothing may be written")
options(optop.revision_suffix = old22)
pass("U22 checkpoint reuse, refusal on identity mismatch, inert without a suffix")

# --- U23: the word-null convention of an input is established, never guessed -------------
tmp23 <- file.path(tempdir(), "mdna_results_MDNA_2015_2016.qs2")   # a REGISTERED name...
obj23 <- list(word = data.table(d_model = 1, d_null = 2))
qs2::qs_save(obj23, tmp23)                                          # ...wrong content
if (!inherits(try(word_null_convention_of(obj23, tmp23), silent = TRUE), "try-error"))
  fail("U23 a file with a legacy NAME but other content must be refused")
tmp23b <- file.path(tempdir(), "some_new_result.qs2"); qs2::qs_save(obj23, tmp23b)
if (!inherits(try(word_null_convention_of(obj23, tmp23b), silent = TRUE), "try-error"))
  fail("U23 an unstamped, unregistered input must be refused as ambiguous")
st23 <- obj23; attr(st23, "run_meta") <- run_meta(list(experiment = "unit"))
if (!identical(word_null_convention_of(st23, tmp23b), "poisson_scored_tokens"))
  fail("U23 a stamped object must answer for itself")
bad23 <- obj23; attr(bad23, "run_meta") <- list(scoring = list(word_null_convention = "??"))
if (!inherits(try(word_null_convention_of(bad23, tmp23b), silent = TRUE), "try-error"))
  fail("U23 an unknown stamp must be refused")
unlink(c(tmp23, tmp23b))
pass("U23 word-null convention: stamp or registered content hash, otherwise an error")

# --- U24: support resolution against a dense evaluation of its definition ---------------
th24 <- lapply(fits, foldin_theta, newdata_dtm = sim_ev15$dtm, seed = 1L)
ph24 <- lapply(fits, phi_from_fit)
part24 <- make_heldout_partition(th24, ph24, sim_ev15$dtm, base15$pi_glob, 1)
k24 <- names(fits)[2L]
r24 <- support_resolution(th24[[k24]], ph24[[k24]], sim_ev15$dtm, part24, as.integer(k24))
N24 <- as.matrix(sim_ev15$dtm); R24 <- part24$rare_mask; P24 <- th24[[k24]] %*% ph24[[k24]]
if (max(abs(r24$observed_pooled_count - rowSums(N24 * R24))) > 1e-10 ||
    max(abs(r24$predicted_pooled_share - rowSums(P24 * R24))) > 1e-12 ||
    !all(r24$retained_cells == rowSums(!R24)))
  fail("U24 support_resolution() disagrees with the dense definition")
s24 <- summarise_resolution(r24)
if (abs(s24$observed_pooled_token_share - sum(N24 * R24) / sum(N24)) > 1e-12 ||
    abs(s24$predicted_pooled_token_share -
          sum(rowSums(N24) * rowSums(P24 * R24)) / sum(N24)) > 1e-12)
  fail("U24 pooled token shares must pool counts, not average document shares")
pass("U24 support resolution == dense definition; pooled shares pool counts")

# --- U25: completion scores the SCORED tokens only (token totals) ------------------------
spl25 <- split_tokens_binomial(sim_ev15$dtm, 0.5, seed = 5L)
if (!isTRUE(all.equal(as.matrix(spl25$foldin + spl25$score), as.matrix(sim_ev15$dtm),
                      check.attributes = FALSE)))
  fail("U25 the token split does not conserve the document counts")
com25 <- suppressMessages(score_heldout(fits, spl25$foldin, spl25$score, base15$pi_glob,
                                        c = 1, metrics = "dev", word_at = 4L))
if (max(abs(com25$doc[K == 4L][match(rownames(spl25$score), doc_id), L] -
            as.numeric(Matrix::rowSums(spl25$score)))) > 0)
  fail("U25 completion document lengths must be those of the scored tokens")
S25 <- as.matrix(spl25$score); B25 <- outer(rowSums(S25), as.numeric(base15$pi_glob))
ok25 <- colSums(S25) > 0 & as.numeric(base15$pi_glob) > 0
ref25 <- 2 * colSums(ifelse(S25 > 0, S25 * log(S25 / B25), 0) - (S25 - B25))
w25 <- com25$word[match(colnames(spl25$score), word_id)]
if (max(abs(w25$d_null[ok25] - ref25[ok25])) > 1e-8 * max(1, max(abs(ref25[ok25]))))
  fail("U25 completion word-level null is not the Poisson form on the scored tokens")
Bfull25 <- outer(rowSums(as.matrix(sim_ev15$dtm)), as.numeric(base15$pi_glob))
wrong25 <- 2 * colSums(ifelse(S25 > 0, S25 * log(S25 / Bfull25), 0) - (S25 - Bfull25))
if (max(abs(w25$d_null[ok25] - wrong25[ok25])) < 1e-6)
  fail("U25 negative control failed: full-document totals give the same null")
pass("U25 completion: scored-token lengths and baseline totals (negative control differs)")

cat("ALL UNIT CHECKS PASSED\n")

