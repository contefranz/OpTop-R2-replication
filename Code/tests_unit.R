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
#   U7 instruments: rows sum to zero; moment test ~ nominal size on iid noise
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
pass("U6 epsilon-rule selects the documented K")

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

# --- U10: Lemma 2 identity (doc-wise == word-wise unbinned fitted deviance) --------
l2 <- lemma2_residual(theta_from_fit(fits[["4"]]), phi_from_fit(fits[["4"]]),
                      sim$dtm)
if (abs(l2$resid) > 1e-8)
  fail("U10 Lemma-2 residual %.3e (dev_doc %.4f vs dev_word %.4f)",
       l2$resid, l2$dev_doc, l2$dev_word)
pass(sprintf("U10 Lemma-2 identity holds (resid %.2e)", l2$resid))

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

cat("ALL UNIT CHECKS PASSED\n")
