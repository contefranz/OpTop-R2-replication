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

#' POINTWISE adjacent rule (exploratory): smallest K whose one-sided upper
#' bound on the adjacent gain is <= eps. Each bound is at level alpha on its
#' own, so scanning the grid carries no familywise guarantee; the certified
#' rules are select_k_total_gain() and select_k_adjacent() below.
#' Returns NA when no grid point qualifies (grid should be extended).
select_k_epsilon <- function(gains_dt, eps = 0.01, alpha = 0.05) {
  eps_v <- eps; alpha_v <- alpha
  stopifnot(all(gains_dt$alpha == alpha_v))
  gains_dt[, {
    # is.na guard: logical indexing with an NA bound would inject NA into `ok`
    # and turn a valid selection at a smaller K into NA
    ok <- K[!is.na(ub_onesided) & ub_onesided <= eps_v]
    .(K_hat = if (length(ok)) min(ok) else NA_integer_,
      eps = eps_v, alpha = alpha_v)
  }, by = metric]
}

# --- Definition 1 (revised): simultaneous bounds and total-gain adequacy ------------

#' Paired per-document gains for a whole comparison family, with simultaneous
#' (Bonferroni) one-sided upper bounds.
#'
#' family = "all_pairs": every K < K' on the grid, M = m(m-1)/2 -- the family
#'   behind total-gain adequacy, A(K) = max_{K' > K} {mu(K') - mu(K)} <= eps;
#' family = "adjacent":  the m - 1 successive steps, M = m - 1 (local adequacy).
#' ub_simul uses z_{1 - alpha/M}; ub_pointwise uses z_{1 - alpha}. The bounds
#' do not involve eps, so several tolerances share one coverage event.
#'
#' The retained-document set must not depend on K (it is defined by D_j(null));
#' this is asserted, because paired gains are meaningless otherwise.
#'
#' @param cluster optional data.frame(doc_id, cluster_id). When given, the SE
#'   of each paired mean is cluster-robust,
#'   sqrt(G/(G-1) * sum_g (sum_{j in g} (d_j - dbar))^2) / n, which keeps the
#'   document-weighted estimand. It addresses dependence WITHIN the evaluation
#'   sample only, not clusters shared between training and evaluation.
paired_gains_all <- function(doc_dt, alpha = 0.05, cluster = NULL,
                             family = c("all_pairs", "adjacent")) {
  family_v <- match.arg(family); alpha_v <- alpha
  stopifnot(all(c("K", "metric", "doc_id", "r2_doc") %in% names(doc_dt)))
  cl_map <- NULL
  if (!is.null(cluster)) {
    cluster <- as.data.table(cluster)
    stopifnot(all(c("doc_id", "cluster_id") %in% names(cluster)),
              !anyDuplicated(cluster$doc_id))
    cl_map <- stats::setNames(as.character(cluster$cluster_id),
                              as.character(cluster$doc_id))
  }
  rows <- list()
  for (m_v in unique(doc_dt$metric)) {
    sub <- doc_dt[metric == m_v, .(doc_id, K, r2_doc)]
    stopifnot(!anyDuplicated(sub, by = c("doc_id", "K")))
    wide <- dcast(sub, doc_id ~ K, value.var = "r2_doc")
    ids <- as.character(wide$doc_id)
    R <- as.matrix(wide[, -1L])
    Ks <- as.integer(colnames(R))
    o <- order(Ks); R <- R[, o, drop = FALSE]; Ks <- Ks[o]
    n_na <- rowSums(is.na(R))
    if (!all(n_na %in% c(0L, ncol(R))))
      stop(sprintf(paste0("metric '%s': the retained-document set varies with K",
                          " (%d documents are NA at some but not all K)"),
                   m_v, sum(!n_na %in% c(0L, ncol(R)))), call. = FALSE)
    keep <- n_na == 0L
    R <- R[keep, , drop = FALSE]; ids <- ids[keep]
    n <- nrow(R); m <- length(Ks)
    if (m < 2L || n < 2L) next
    g <- NULL
    if (!is.null(cl_map)) {
      g <- cl_map[ids]
      if (anyNA(g)) stop("cluster_id missing for some retained documents",
                         call. = FALSE)
    }
    pr <- if (family_v == "all_pairs") t(utils::combn(m, 2L)) else
      cbind(seq_len(m - 1L), seq_len(m - 1L) + 1L)
    M <- nrow(pr)
    mu <- sdv <- se <- numeric(M)
    for (i in seq_len(M)) {
      v <- R[, pr[i, 2L]] - R[, pr[i, 1L]]
      mu[i] <- mean(v); sdv[i] <- stats::sd(v)
      se[i] <- if (is.null(g)) sdv[i] / sqrt(n) else {
        s <- rowsum(v - mu[i], g)[, 1L]
        G <- length(s)
        sqrt(G / (G - 1) * sum(s^2)) / n
      }
    }
    z_sim <- stats::qnorm(1 - alpha_v / M); z_pt <- stats::qnorm(1 - alpha_v)
    rows[[m_v]] <- data.table(
      metric = m_v, K = Ks[pr[, 1L]], K_to = Ks[pr[, 2L]],
      delta_mean = mu, delta_sd = sdv, se = se, n = n,
      n_clusters = if (is.null(g)) NA_integer_ else length(unique(g)),
      M = M, z_crit = z_sim,
      ub_simul = mu + z_sim * se, ub_pointwise = mu + z_pt * se,
      alpha = alpha_v, se_type = if (is.null(g)) "iid" else "cluster",
      family = family_v)
  }
  rbindlist(rows)
}

#' Shared selector core. A candidate is certified only when EVERY bound behind
#' it is finite and <= eps; the grid maximum is never a candidate (it has no
#' larger model to be compared with). No certified candidate -> K_hat = NA and
#' certified = FALSE: the rule reports "none certified", never the boundary.
.select_certified <- function(per_k, eps_v, rule_v, meta) {
  per_k[, ok := !is.na(max_ub) & is.finite(max_ub) & max_ub <= eps_v]
  out <- per_k[, {
    i <- which(ok)
    if (length(i)) {
      j <- i[which.min(K[i])]
      .(K_hat = as.integer(K[j]), certified = TRUE,
        n_comparisons = as.integer(n_cmp[j]), max_ub = max_ub[j])
    } else {
      .(K_hat = NA_integer_, certified = FALSE,
        n_comparisons = NA_integer_, max_ub = NA_real_)
    }
  }, by = metric]
  out <- merge(out, meta, by = "metric", sort = FALSE)
  out[, `:=`(eps = eps_v, rule = rule_v)]
  out[]
}

#' Total-gain epsilon-adequacy (primary rule): smallest non-boundary K such that
#' the simultaneous upper bound on the gain to EVERY larger candidate is <= eps.
#' `n_comparisons` is the number of larger candidates behind the certificate; a
#' selection at the penultimate grid point rests on a single comparison.
select_k_total_gain <- function(pairs_dt, eps = 0.01) {
  stopifnot(nrow(pairs_dt) > 0L, all(pairs_dt$family == "all_pairs"))
  per_k <- pairs_dt[, .(max_ub = if (all(is.finite(ub_simul))) max(ub_simul)
                                 else NA_real_,
                        n_cmp = .N), by = .(metric, K)]
  meta <- unique(pairs_dt[, .(metric, alpha, M, se_type)])
  .select_certified(per_k, eps, "total_gain", meta)
}

#' Adjacent-gain adequacy from paired_gains_all(family = "adjacent"):
#' simultaneous = TRUE  -> local rule with Bonferroni over the m - 1 steps;
#' simultaneous = FALSE -> the original pointwise rule (exploratory).
#' Either certifies the NEXT grid step only, never the total remaining gain.
select_k_adjacent <- function(adj_dt, eps = 0.01, simultaneous = TRUE) {
  stopifnot(nrow(adj_dt) > 0L, all(adj_dt$family == "adjacent"))
  ub <- if (simultaneous) adj_dt$ub_simul else adj_dt$ub_pointwise
  per_k <- data.table(metric = adj_dt$metric, K = adj_dt$K, max_ub = ub,
                      n_cmp = 1L)
  meta <- unique(adj_dt[, .(metric, alpha, M, se_type)])
  if (!simultaneous) meta[, M := 1L]
  .select_certified(per_k, eps,
                    if (simultaneous) "adjacent_simultaneous"
                    else "adjacent_pointwise", meta)
}

#' Simultaneous adjacent rule from a CACHED paired_gains() table (columns
#' metric, K, delta_mean, se), so results saved without document-level scores
#' (E1, E2, the second DGP) can be re-selected without rescoring.
select_k_adjacent_simultaneous <- function(gains_dt, eps = 0.01, alpha = 0.05) {
  alpha_v <- alpha
  stopifnot(all(c("metric", "K", "delta_mean", "se") %in% names(gains_dt)),
            !anyDuplicated(gains_dt, by = c("metric", "K")))
  g <- copy(gains_dt)[, M := .N, by = metric]
  per_k <- g[, .(metric, K, n_cmp = 1L,
                 max_ub = delta_mean + stats::qnorm(1 - alpha_v / M) * se)]
  meta <- unique(g[, .(metric, M)])[, `:=`(alpha = alpha_v, se_type = "iid")]
  .select_certified(per_k, eps, "adjacent_simultaneous", meta)
}

#' All three rules at several tolerances from one tidy document table.
#' A metric with fewer than two retained documents (e.g. a support that has
#' collapsed, so every null discrepancy sits below the floor) cannot certify
#' anything: it is reported as "none certified" with n_retained = 0 rather
#' than dropped or raised as an error.
select_k_all_rules <- function(doc_dt, eps_grid = c(0.01, 0.005), alpha = 0.05,
                               cluster = NULL) {
  RULES <- c("adjacent_pointwise", "adjacent_simultaneous", "total_gain")
  adj <- paired_gains_all(doc_dt, alpha, cluster, "adjacent")
  allp <- paired_gains_all(doc_dt, alpha, cluster, "all_pairs")
  res <- if (nrow(adj)) rbindlist(lapply(eps_grid, function(e) rbindlist(list(
    select_k_adjacent(adj, e, simultaneous = FALSE),
    select_k_adjacent(adj, e, simultaneous = TRUE),
    select_k_total_gain(allp, e)), use.names = TRUE))) else NULL
  if (!is.null(res))
    res <- merge(res, unique(adj[, .(metric, n_retained = n)]), by = "metric",
                 sort = FALSE)
  lost <- setdiff(unique(doc_dt$metric), unique(adj$metric))
  if (length(lost)) {
    none <- CJ(metric = lost, eps = eps_grid, rule = RULES, sorted = FALSE)[
      , `:=`(K_hat = NA_integer_, certified = FALSE, n_comparisons = NA_integer_,
             max_ub = NA_real_, alpha = alpha, M = NA_integer_,
             se_type = NA_character_, n_retained = 0L)]
    res <- rbindlist(list(res, none), use.names = TRUE)
  }
  # fixed layout: callers stack these tables, and "none certified" rows are
  # built separately from certified ones
  setcolorder(res, c("metric", "rule", "eps", "K_hat", "certified",
                     "n_comparisons", "max_ub", "M", "alpha", "se_type",
                     "n_retained"))
  res[]
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
#'
#' `word_d_model` (optional): the PACKAGE's word-level fitted deviances for ALL
#' W words, in vocabulary order (score_insample(..., word_at = K)$word$d_model).
#' Without it the function can only re-add one vector of cell terms in two
#' orders, which agrees to rounding by construction and verifies nothing about
#' the implementation (the pre-revision check did exactly that). With it, the
#' directly computed document-wise total is compared with a separately
#' implemented code path, in total and word by word:
#'   resid_pkg     = dev_doc - sum_w D_w^pkg
#'   max_word_diff = max_w | D_w^pkg - D_w^direct |,
#' with D_w^direct the Poisson-corrected word deviance of Section 2.4,
#'   D_w = 2 sum_j [ N_jw log(N_jw / E_jw) - (N_jw - E_jw) ].
#'
#' ZERO-PROBABILITY CELLS. The lemma assumes E_jw > 0 wherever N_jw > 0. An
#' unsmoothed topic-word matrix (WarpLDA returns count-normalised phi with most
#' entries exactly zero) can leave a few observed cells with fitted probability
#' exactly 0, where the unbinned deviance is infinite. The package floors the
#' EXPECTED COUNT at 1e-12 on its word-level path; `e_floor` reproduces that
#' convention so the two code paths are compared like for like, and the number
#' and deviance share of floored cells are returned because every word-level
#' quantity touching such a cell depends on this arbitrary constant. (The
#' document-level indices are immune: a zero-probability word always falls in
#' the min-bin of the harmonised support.)
lemma2_residual <- function(theta, phi, dtm, word_d_model = NULL,
                            e_floor = 1e-12) {
  L <- as.numeric(Matrix::rowSums(dtm))
  P <- theta %*% phi                             # J x W fitted probabilities
  Tm <- as(as(dtm, "generalMatrix"), "TsparseMatrix")
  ii <- Tm@i + 1L; jj <- Tm@j + 1L
  E_raw <- P[cbind(ii, jj)] * L[ii]
  floored <- E_raw < e_floor
  E_nz <- pmax(E_raw, e_floor)
  cell <- 2 * Tm@x * log(Tm@x / E_nz)            # 2 N log(N/E), N > 0 only
  by_doc  <- rowsum(cell, ii)[, 1L]              # genuine per-document sums
  by_word <- rowsum(cell, jj)                    # genuine per-word sums
  dev_doc <- sum(by_doc); dev_word <- sum(by_word[, 1L])
  out <- list(dev_doc = dev_doc, dev_word = dev_word,
              resid = dev_doc - dev_word,
              n_cells = length(cell), n_floored = sum(floored),
              token_share_floored = sum(Tm@x[floored]) / sum(Tm@x),
              dev_share_floored = sum(cell[floored]) / dev_doc,
              n_words_floored = length(unique(jj[floored])))
  if (!is.null(word_d_model)) {
    stopifnot(length(word_d_model) == ncol(dtm))
    dw <- numeric(ncol(dtm))
    dw[as.integer(rownames(by_word))] <- by_word[, 1L]
    E_w <- as.numeric(crossprod(P, L))           # sum_j L_j p_jw
    N_w <- as.numeric(Matrix::colSums(dtm))
    dw_direct <- dw - 2 * (N_w - E_w)            # Poisson correction
    out$dev_word_pkg <- sum(word_d_model)
    out$resid_pkg <- dev_doc - out$dev_word_pkg
    out$rel_resid_pkg <- out$resid_pkg / dev_doc
    out$max_word_diff <- max(abs(word_d_model - dw_direct))
  }
  out
}
