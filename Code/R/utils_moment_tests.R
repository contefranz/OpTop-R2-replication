# =============================================================================
# utils_moment_tests.R
# Section 4: held-out moment-based specification tests.
#
#   Test 1  frequency-contrast screen   (q = 1;  t -> N(0,1), W1 -> chi2_1)
#   Test 2  frequency-stratified Wald   (q = B-1; chi2_{B-1})
#   Test 3  fit-stratified Wald         (q = S-1; chi2_{S-1})
#
# Instruments are built from the TRAINING sample only and each row sums to
# zero over the vocabulary (eq. 53). Held-out residuals are
# e_j = d_j - i_hat_j^{K,tr} on the training-support vocabulary (probability
# scale), so moments are in probability-mass units (the paper's effect size).
# =============================================================================

library(data.table)
library(Matrix)

# --- Instruments ----------------------------------------------------------------

#' Frequency strata from training corpus frequency (deterministic ranking).
#' Returns an integer stratum id per word, 1 = lowest frequency, B = highest.
.freq_strata <- function(train_dtm, B = 5L) {
  f <- as.numeric(Matrix::colSums(train_dtm))
  as.integer(cut(rank(f, ties.method = "first"),
                 breaks = quantile(rank(f, ties.method = "first"),
                                   probs = seq(0, 1, length.out = B + 1L)),
                 include.lowest = TRUE, labels = FALSE))
}

#' Test 1 instrument: top vs bottom frequency quintile (1 x W, sums to zero).
make_instrument_freq_contrast <- function(train_dtm, n_strata = 5L) {
  s <- .freq_strata(train_dtm, n_strata)
  hi <- s == n_strata; lo <- s == 1L
  z <- numeric(ncol(train_dtm))
  z[hi] <- 1 / sum(hi); z[lo] <- -1 / sum(lo)
  Z <- matrix(z, nrow = 1L,
              dimnames = list("hi_vs_lo", colnames(train_dtm)))
  Z
}

#' Test 2 instruments: stratum b vs highest-frequency reference ((B-1) x W).
make_instruments_freq_strata <- function(train_dtm, B = 5L) {
  s <- .freq_strata(train_dtm, B)
  W <- ncol(train_dtm)
  Z <- matrix(0, nrow = B - 1L, ncol = W,
              dimnames = list(paste0("f", seq_len(B - 1L), "_vs_f", B),
                              colnames(train_dtm)))
  nB <- sum(s == B)
  for (b in seq_len(B - 1L)) {
    Z[b, s == b] <- 1 / sum(s == b)
    Z[b, s == B] <- -1 / nB
  }
  Z
}

#' Test 3 instruments: strata of a TRAINING word-level fit score (S-1) x W;
#' words failing the filter get zero entries (rows still sum to zero).
#' Reference stratum = highest-score (best-fit) words.
make_instruments_fit_strata <- function(word_scores, vocab, S = 5L,
                                        min_docfreq = 5L, docfreq = NULL) {
  sc <- word_scores[match(vocab, word_id)]
  keep <- !is.na(sc$r2_word)
  if (!is.null(docfreq)) keep <- keep & docfreq >= min_docfreq
  W <- length(vocab)
  strat <- rep(NA_integer_, W)
  strat[keep] <- as.integer(cut(rank(sc$r2_word[keep], ties.method = "first"),
                                breaks = quantile(rank(sc$r2_word[keep],
                                                       ties.method = "first"),
                                                  probs = seq(0, 1, length.out = S + 1L)),
                                include.lowest = TRUE, labels = FALSE))
  Z <- matrix(0, nrow = S - 1L, ncol = W,
              dimnames = list(paste0("s", seq_len(S - 1L), "_vs_s", S), vocab))
  nS <- sum(strat == S, na.rm = TRUE)
  for (b in seq_len(S - 1L)) {
    Z[b, which(strat == b)] <- 1 / sum(strat == b, na.rm = TRUE)
    Z[b, which(strat == S)] <- -1 / nS
  }
  Z
}

# --- Residuals -------------------------------------------------------------------

#' Held-out residual matrix e_j = d_j - theta_j' Phi on the training vocabulary.
resid_heldout <- function(theta, phi, score_dtm) {
  L <- Matrix::rowSums(score_dtm)
  D <- as.matrix(score_dtm / L)
  D - theta %*% phi
}

# --- Wald machinery -----------------------------------------------------------------

#' Wald test of H0: E[Z e_j] = 0 from a J x q moment matrix G = E Z'.
#' Guards the covariance inversion (condition number -> pseudo-inverse + flag).
#'
#' `center` (optional q-vector): ALSO test H0: mu = center -- the conditional
#' truth under the training fit (estimated from a large evaluation set). Per
#' Remark 8, the raw mu = 0 null is rejected at scale by ANY residual
#' imbalance of the fitted model (e.g. VEM smoothing bias, ~1e-7 mass units);
#' the centered statistic isolates the calibration of the Wald machinery.
moment_test <- function(Z, E, test_label = "test", center = NULL) {
  moment_test_from_G(E %*% t(Z), test_label, center)
}

#' Wald machinery on a precomputed J x q moment matrix G. Split out so large
#' evaluation sets can accumulate G block-wise (G_block = E_block %*% t(Z))
#' without materializing the dense J x W residual matrix (run_mdna.R).
moment_test_from_G <- function(G, test_label = "test", center = NULL) {
  G <- as.matrix(G)
  J <- nrow(G); q <- ncol(G)
  gbar <- colMeans(G)
  S <- stats::cov(G)
  flag <- ""
  Sinv <- tryCatch({
    if (q == 1L) {
      matrix(1 / S, 1L, 1L)
    } else if (rcond(S) < 1e-12) {
      flag <- "pseudo-inverse"
      MASS::ginv(S)
    } else {
      solve(S)
    }
  }, error = function(e) { flag <<- "pseudo-inverse"; MASS::ginv(S) })
  wald <- function(v) as.numeric(J * t(v) %*% Sinv %*% v)
  W <- wald(gbar)
  res <- data.table(
    test = test_label, stat = W, df = q,
    pval = pchisq(W, df = q, lower.tail = FALSE),
    gbar_absmax = max(abs(gbar)), J_ev = J, flag = flag
  )
  if (!is.null(center)) {
    Wc <- wald(gbar - center)
    res[, `:=`(stat_centered = Wc,
               pval_centered = pchisq(Wc, df = q, lower.tail = FALSE))]
  }
  t_strata <- sqrt(J) * gbar / sqrt(diag(S))
  p_strata <- 2 * pnorm(-abs(t_strata))
  list(
    result = res,
    strata = data.table(
      test = test_label, stratum = colnames(G) %||% paste0("m", seq_len(q)),
      gbar = gbar, t = t_strata, pval = p_strata,
      pval_bh = p.adjust(p_strata, method = "BH")
    ),
    gbar = gbar
  )
}

#' Build the instrument set once per (training corpus, K).
make_instrument_set <- function(train_dtm, word_scores_train = NULL,
                                B = 5L, S = 5L, min_docfreq = 5L) {
  Zs <- list(
    T1_freq_contrast = make_instrument_freq_contrast(train_dtm),
    T2_freq_strata   = make_instruments_freq_strata(train_dtm, B = B)
  )
  if (!is.null(word_scores_train)) {
    docfreq_tr <- as.integer(Matrix::colSums(train_dtm > 0))
    Zs$T3_fit_strata <- make_instruments_fit_strata(
      word_scores_train, colnames(train_dtm), S = S,
      min_docfreq = min_docfreq, docfreq = docfreq_tr)
  }
  Zs
}

#' Run the Section-4 battery on a residual matrix (theta already folded in).
#' `centers`: optional named list (by test) of conditional-truth vectors.
run_moment_battery <- function(Zs, E, centers = NULL) {
  out <- lapply(names(Zs), function(nm)
    moment_test(Zs[[nm]], E, nm, center = centers[[nm]]))
  names(out) <- names(Zs)
  list(results = rbindlist(lapply(out, `[[`, "result"), fill = TRUE),
       strata = rbindlist(lapply(out, `[[`, "strata")),
       gbars = lapply(out, `[[`, "gbar"))
}

#' Convenience wrapper: fold in theta, then run the battery (single use;
#' drivers that also score R2 should fold in once and call run_moment_battery).
run_moment_tests <- function(fit, foldin_dtm, score_dtm, train_dtm,
                             word_scores_train = NULL,
                             B = 5L, S = 5L, min_docfreq = 5L,
                             foldin_seed = 1L, centers = NULL) {
  theta <- foldin_theta(fit, foldin_dtm, seed = foldin_seed)
  phi <- phi_from_fit(fit)
  E <- resid_heldout(theta, phi, score_dtm)
  Zs <- make_instrument_set(train_dtm, word_scores_train, B, S, min_docfreq)
  run_moment_battery(Zs, E, centers)
}
