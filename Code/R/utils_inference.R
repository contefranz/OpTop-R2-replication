# =============================================================================
# utils_inference.R
# Cross-document inference layer of Section 3:
#   * Proposition 2:  Macro CI, paired adjacent-gain tests;
#   * Definition 1:   epsilon-adequacy selection rule K_hat_{D,eps,alpha};
#   * Remark 6:       delta-method SE for the Micro-Macro gap;
#   * Proposition 1iii: exact length / atypicality / interaction decomposition;
#   * Section 3.8:    word-level filtering rule (made precise per referee).
#
# All functions consume the tidy per-document tables produced by
# score_heldout()/score_insample(): columns K, metric, doc_id, r2_doc,
# d_model, d_null, L. NA r2_doc marks degenerate documents (D_j(null) = 0),
# which are excluded per the paper's J_{ev,+} convention.
# =============================================================================

library(data.table)

# --- Proposition 2: Macro index CI --------------------------------------------------

macro_ci <- function(doc_dt, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  out <- doc_dt[!is.na(r2_doc), .(
    r2_macro = mean(r2_doc),
    sd = sd(r2_doc),
    n = .N
  ), by = .(K, metric)]
  out[, se := sd / sqrt(n)]
  out[, `:=`(lwr = r2_macro - z * se, upr = r2_macro + z * se, level = level)]
  out[]
}

# --- Proposition 2 (iv): paired adjacent gains --------------------------------------

#' Paired per-document differences between adjacent grid points.
#' The gain indexed at K is Delta(K) = R2_j(succ(K)) - R2_j(K); the one-sided
#' upper bound at level alpha feeds the Definition-1 rule.
paired_gains <- function(doc_dt, alpha = 0.05, level = 0.95) {
  alpha_v <- alpha
  z1 <- qnorm(1 - alpha_v)
  z2 <- qnorm(1 - (1 - level) / 2)
  Ks <- sort(unique(doc_dt$K))
  if (length(Ks) < 2L) return(data.table())
  rows <- vector("list", length(Ks) - 1L)
  for (i in seq_len(length(Ks) - 1L)) {
    a <- doc_dt[K == Ks[i], .(metric, doc_id, r2_a = r2_doc)]
    b <- doc_dt[K == Ks[i + 1L], .(metric, doc_id, r2_b = r2_doc)]
    ab <- merge(a, b, by = c("metric", "doc_id"))[!is.na(r2_a) & !is.na(r2_b)]
    rows[[i]] <- ab[, .(
      K = Ks[i], K_next = Ks[i + 1L],
      delta_mean = mean(r2_b - r2_a),
      delta_sd = sd(r2_b - r2_a),
      n = .N
    ), by = metric]
  }
  out <- rbindlist(rows)
  out[, se := delta_sd / sqrt(n)]
  out[, `:=`(
    ub_onesided = delta_mean + z1 * se,
    lwr = delta_mean - z2 * se,
    upr = delta_mean + z2 * se,
    alpha = alpha_v
  )]
  out[]
}

# --- Definition 1: epsilon-adequacy selection ---------------------------------------

#' Smallest K whose one-sided upper bound on the adjacent gain is <= eps.
#' Returns NA when no grid point qualifies (grid should be extended).
select_k_epsilon <- function(gains_dt, eps = 0.01, alpha = 0.05) {
  eps_v <- eps; alpha_v <- alpha
  stopifnot(all(gains_dt$alpha == alpha_v))
  gains_dt[, {
    ok <- K[ub_onesided <= eps_v]
    .(K_hat = if (length(ok)) min(ok) else NA_integer_,
      eps = eps_v, alpha = alpha_v)
  }, by = metric]
}

# --- Remark 6: Micro-Macro gap with delta-method SE ---------------------------------

#' gap = Cov_hat(r_j, D_j(null)) / mean(D_j(null)) over J_+ (exact identity,
#' Prop 1ii). SE from the influence function of the smooth functional.
micro_macro_gap_ci <- function(doc_dt, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  out <- doc_dt[!is.na(r2_doc) & d_null > 0, {
    u <- r2_doc; v <- d_null; n <- .N
    ub <- mean(u); vb <- mean(v)
    suv <- mean((u - ub) * (v - vb))
    gap <- suv / vb
    psi <- ((u - ub) * (v - vb) - suv) / vb - gap * (v - vb) / vb
    .(gap = gap, se = sd(psi) / sqrt(n), n = n)
  }, by = .(K, metric)]
  out[, `:=`(lwr = gap - z * se, upr = gap + z * se, level = level)]
  out[]
}

# --- Proposition 1(iii): exact gap decomposition (Deviance family) -------------------

#' D_j(null) = 2 L_j kappa_j, so kappa_j = d_null / (2 L_j). The gap decomposes
#' exactly into a length channel, an atypicality channel, and an interaction
#' channel (eq. 31); `resid` checks the identity at machine precision.
gap_decomposition <- function(doc_dt_dev) {
  stopifnot(all(doc_dt_dev$metric == "dev"))
  out <- doc_dt_dev[!is.na(r2_doc) & d_null > 0, {
    r <- r2_doc; L_ <- L; kap <- d_null / (2 * L_)
    n <- .N
    cv <- function(a, b) mean((a - mean(a)) * (b - mean(b)))
    Lk <- L_ * kap
    denom <- mean(Lk)
    gap <- cv(r, Lk) / denom
    ch_len <- mean(kap) * cv(r, L_) / denom
    ch_aty <- mean(L_) * cv(r, kap) / denom
    ch_int <- cv(r, (L_ - mean(L_)) * (kap - mean(kap))) / denom
    .(gap = gap, ch_length = ch_len, ch_atypicality = ch_aty,
      ch_interaction = ch_int, resid = gap - (ch_len + ch_aty + ch_int), n = n)
  }, by = K]
  out[]
}

# --- Section 3.8: word-level filtering rule ------------------------------------------

#' Precise rule (referee Minor 5): keep word w iff
#'   (i)  document frequency in the SCORED corpus >= min_docfreq, and
#'   (ii) expected corpus count under the null baseline >= min_expected
#'        (B_w = pi_glob(w) * sum_j L_j), and
#'   (iii) D_w(null) > 0 (non-degenerate).
word_index_filtered <- function(word_dt, score_dtm, pi_glob,
                                min_docfreq = 5L, min_expected = 5) {
  docfreq <- Matrix::colSums(score_dtm > 0)
  B_w <- as.numeric(pi_glob) * sum(Matrix::rowSums(score_dtm))
  meta <- data.table(word_id = colnames(score_dtm),
                     doc_freq = as.integer(docfreq), B_w = B_w)
  out <- merge(word_dt, meta, by = "word_id", all.x = TRUE, sort = FALSE)
  out[, keep := doc_freq >= min_docfreq & B_w >= min_expected & d_null > 0]
  summary <- out[keep == TRUE, .(
    r2_micro_word = 1 - sum(d_model) / sum(d_null),
    r2_macro_word = mean(r2_word, na.rm = TRUE),
    n_words = .N
  ), by = .(K, metric)]
  list(word = out[], summary = summary,
       rule = list(min_docfreq = min_docfreq, min_expected = min_expected))
}

#' Word-level Micro/Macro indices and their divergence across the K grid
#' (Section 3.8.2-3.8.3). The gap w-Micro - w-Macro is the frequency-space
#' analogue of the document-level Micro-Macro length-bias gap (Section 3.6):
#' a large positive gap means the model fits common (high-baseline-discrepancy)
#' words well but much of the rarer vocabulary poorly.
#'
#' @param word_dt tidy word table (K, metric, word_id, r2_word, d_model, d_null)
#'   for one evaluation mode.
word_micro_macro_over_k <- function(word_dt, score_dtm, pi_glob,
                                    min_docfreq = 5L, min_expected = 5) {
  wf <- word_index_filtered(word_dt, score_dtm, pi_glob,
                            min_docfreq, min_expected)
  out <- wf$summary[, .(K, metric, r2_micro_word, r2_macro_word,
                        gap = r2_micro_word - r2_macro_word, n_words)]
  list(curve = out[order(metric, K)], word = wf$word, rule = wf$rule)
}

#' Lemma 2 identity (eq. 51): on the UNBINNED vocabulary the total fitted
#' deviance summed over documents equals the total summed over words,
#' 2 sum_jw N_jw log(N_jw / E_jw). Returns doc-wise total, word-wise total,
#' and their residual (must be ~0). E = theta %*% phi scaled by L_j; the
#' identity holds cell-wise, so both aggregations are computed from the same
#' sparse N and dense E without any rare-word binning.
lemma2_residual <- function(theta, phi, dtm) {
  L <- as.numeric(Matrix::rowSums(dtm))
  E <- (theta %*% phi) * L                       # J x W expected counts
  Tm <- as(as(dtm, "generalMatrix"), "TsparseMatrix")
  idx <- cbind(Tm@i + 1L, Tm@j + 1L)
  cell <- 2 * Tm@x * log(Tm@x / pmax(E[idx], 1e-300))   # N log(N/E), N>0 only
  dev_doc  <- sum(cell[order(Tm@i)])             # aggregate by document
  dev_word <- sum(cell[order(Tm@j)])             # aggregate by word
  list(dev_doc = dev_doc, dev_word = dev_word,
       resid = dev_doc - dev_word)
}
