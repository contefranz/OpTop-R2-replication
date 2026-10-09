# =============================================================================
# utils_heldout.R
# Held-out evaluation exactly as in Section 3.7 of the paper, plus the native
# in-sample scoring path.
#
# Held-out protocol:
#   * global objects (Phi_hat, pi_glob) come from the TRAINING sample only;
#   * theta_hat for evaluation documents is folded in from fold-in tokens
#     (reconstruction: all tokens; completion: the fold-in half only);
#   * the harmonized rare set uses the training baseline:
#       rare[j, w]  <=>  min( pi_glob^tr(w), min_K i^{K,tr}_jw ) < c / L_j,
#     with L_j the length of the SCORED tokens;
#   * discrepancies are evaluated by OpTop's exported index functions, reached
#     through the official `nlp_topic_fit` adapter (see as_pseudo_fit()).
#
# Since OpTop 0.14.0, optop_make_partition() accepts an EXTERNAL baseline via
# `pi_glob=` and computes this exact rule in compiled (OpenMP-blocked) cores, so
# the production partition is delegated to the package; the pure-R reference
# implementation is kept below as .make_heldout_partition_ref() and equality is
# enforced by U15 in the test suite (plus the in-sample crosscheck in U1).
# One report field is recomputed locally: the package defines
# chisq_min_report$excluded_mass as a per-document mean, while this pipeline's
# T7b/min-bin tables use the pooled corpus-share definition.
# =============================================================================

library(OpTop)
library(Matrix)
library(data.table)

# --- Held-out harmonized partition (Section 3.7 rule) -----------------------------

#' Production builder: native OpTop core with the TRAINING baseline passed in.
#' Same signature and return shape as the pure-R reference below.
make_heldout_partition <- function(theta_list, phi_list, score_dtm, pi_glob,
                                   c = 1) {
  stopifnot(length(theta_list) == length(phi_list), length(theta_list) >= 1L)
  pfs <- lapply(seq_along(theta_list), function(k)
    as_pseudo_fit(theta_list[[k]], phi_list[[k]],
                  doc_ids = rownames(score_dtm), vocab = colnames(score_dtm)))
  part <- OpTop::optop_make_partition(models = pfs, dtm = score_dtm, c = c,
                                      n_threads = 1L,   # corpus-level parallelism upstream
                                      pi_glob = as.numeric(pi_glob))
  part <- .densify_partition(part, score_dtm)
  dimnames(part$rare_mask) <- dimnames(score_dtm)

  # pooled corpus-share definition of excluded_mass (pipeline convention;
  # the package reports a per-document mean instead)
  has_min <- rowSums(part$rare_mask) > 0
  excluded <- has_min & !part$chisq_min_ok
  Tm <- as(as(score_dtm, "generalMatrix"), "TsparseMatrix")
  hit <- part$rare_mask[cbind(Tm@i + 1L, Tm@j + 1L)]
  tot <- sum(Tm@x)
  part$chisq_min_report$excluded_mass <-
    if (tot > 0) sum(Tm@x[hit & excluded[Tm@i + 1L]]) / tot else 0
  part
}

#' OpTop >= 0.20 returns the partition in compressed form ("format 2"):
#' per-document 0-based column indices of the NON-rare words, with no dense
#' rare_mask. This helper rebuilds the dense logical mask (TRUE = rare) that
#' the pipeline consumes downstream; a toy-verified no-op for <= 0.14 objects.
#' Numerical equality with the frozen pure-R reference stays gated by U15.
.densify_partition <- function(part, score_dtm) {
  if (!is.null(part$rare_mask)) return(part)
  stopifnot(!is.null(part$nonrare_offsets), !is.null(part$nonrare_words),
            identical(as.character(part$vocab), colnames(score_dtm)))
  J <- nrow(score_dtm); W <- ncol(score_dtm)
  n_per_doc <- diff(part$nonrare_offsets)
  stopifnot(length(n_per_doc) == J,
            length(part$nonrare_words) == sum(n_per_doc))
  rm_ <- matrix(TRUE, J, W)
  rm_[cbind(rep.int(seq_len(J), n_per_doc), part$nonrare_words + 1L)] <- FALSE
  part$rare_mask <- rm_
  part
}

#' Pure-R reference implementation (Section 3.7 rule, pre-0.14.0 production
#' path). Retained ONLY for the U15 native==reference gate; not called by the
#' pipeline.
.make_heldout_partition_ref <- function(theta_list, phi_list, score_dtm,
                                        pi_glob, c = 1) {
  stopifnot(length(theta_list) == length(phi_list), length(theta_list) >= 1L)
  J <- nrow(score_dtm); W <- ncol(score_dtm)
  L <- as.numeric(Matrix::rowSums(score_dtm))
  tau <- c / L

  # pass 1: elementwise min over the grid of fitted word probabilities
  M <- NULL
  for (k in seq_along(theta_list)) {
    I_k <- theta_list[[k]] %*% phi_list[[k]]
    M <- if (is.null(M)) I_k else pmin(M, I_k)
  }

  # rare  <=>  min(pi_glob, M) < tau  (tau recycles down columns: J-vector)
  rare_mask <- (M < tau) | outer(tau, as.numeric(pi_glob), FUN = ">")
  dimnames(rare_mask) <- dimnames(score_dtm)

  # pass 2: min-bin masses for the Pearson inclusion rule (as in the package)
  rare_mass_min <- NULL
  for (k in seq_along(theta_list)) {
    I_k <- theta_list[[k]] %*% phi_list[[k]]
    m_k <- rowSums(I_k * rare_mask)
    rare_mass_min <- if (is.null(rare_mass_min)) m_k else pmin(rare_mass_min, m_k)
  }
  E_min_min <- L * rare_mass_min
  B_min <- L * as.numeric(rare_mask %*% as.numeric(pi_glob))
  chisq_min_ok <- pmin(E_min_min, B_min) >= c

  has_min <- rowSums(rare_mask) > 0
  excluded <- has_min & !chisq_min_ok
  Tm <- as(as(score_dtm, "generalMatrix"), "TsparseMatrix")
  hit <- rare_mask[cbind(Tm@i + 1L, Tm@j + 1L)]
  tot <- sum(Tm@x)
  excluded_mass <- if (tot > 0) sum(Tm@x[hit & excluded[Tm@i + 1L]]) / tot else 0

  list(
    rare_mask = rare_mask,
    L = L,
    chisq_min_ok = chisq_min_ok,
    chisq_min_report = list(
      n_excluded = sum(excluded),
      share = mean(excluded),
      excluded_mass = excluded_mass
    ),
    c = c
  )
}

# --- Scoring: shared assembly ------------------------------------------------------

.metric_fun <- function(metric) {
  switch(metric,
    dev   = OpTop::optop_index_deviance,
    chisq = OpTop::optop_index_chisq,
    se    = OpTop::optop_index_se,
    stop("unknown metric: ", metric)
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.collect_doc <- function(res, K, metric, L, ids) {
  data.table(
    K = K, metric = metric,
    doc_id = as.character(names(res$r2_doc) %||% ids),
    r2_doc = as.numeric(res$r2_doc),
    d_model = as.numeric(res$d_model),
    d_null = as.numeric(res$d_null),
    L = L
  )
}

#' `d_null_override` (held-out Deviance only): the baseline word deviance in the
#' Poisson form of Section 2.4, from .word_null_dev_poisson(); r2_word is then
#' recomputed from it.
.collect_word <- function(res, K, metric, ids, d_null_override = NULL) {
  d_model <- as.numeric(res$d_model)
  d_null <- as.numeric(res$d_null); r2 <- as.numeric(res$r2_word)
  if (!is.null(d_null_override)) {
    d_null <- as.numeric(d_null_override)
    r2 <- ifelse(d_null > 0, 1 - d_model / d_null, NA_real_)
  }
  data.table(
    K = K, metric = metric,
    word_id = as.character(names(res$r2_word) %||% ids),
    r2_word = r2, d_model = d_model, d_null = d_null
  )
}

#' Baseline word-level deviance in POISSON form,
#'   D_w(null) = 2 sum_j [ N_jw log(N_jw / B_jw) - (N_jw - B_jw) ],  B_jw = L_j pi_w,
#' as Section 2.4 of the paper defines it (every summand is nonnegative).
#'
#' Why it is computed here and not taken from the package: OpTop <= 0.20.1 (0.20.1 adds a validity warning only)
#' returns the word-level FITTED deviance with its linear term but the word-level
#' NULL deviance without it (2 sum_j N log(N/B); src/index_core.cpp). In-sample
#' the omitted term is identically zero, because the baseline is the marginal of
#' the scored corpus (sum_j B_jw = N_w), so nothing in-sample -- including the
#' training word scores behind Test 3 -- is affected. HELD-OUT, the baseline
#' comes from training and the counts from evaluation, so N_w - B_w != 0: the
#' term carries the contribution of the documents in which a word does NOT
#' occur, and dropping it makes words concentrated in few documents look far
#' worse fitted than the baseline. Computing the null directly makes the
#' pipeline correct under either package version. Expected counts are floored
#' at `eps` inside the logarithm only, as in the package.
.word_null_dev_poisson <- function(score_dtm, pi_glob, eps = 1e-12) {
  L <- as.numeric(Matrix::rowSums(score_dtm)); pi_w <- as.numeric(pi_glob)
  stopifnot(length(pi_w) == ncol(score_dtm))
  Tm <- as(as(score_dtm, "generalMatrix"), "TsparseMatrix")
  jj <- Tm@j + 1L
  b <- pmax(pi_w[jj] * L[Tm@i + 1L], eps)
  d <- numeric(ncol(score_dtm))
  if (length(jj)) {
    agg <- rowsum(2 * Tm@x * log(Tm@x / b), jj)
    d[as.integer(rownames(agg))] <- agg[, 1L]
  }
  d - 2 * (as.numeric(Matrix::colSums(score_dtm)) - pi_w * sum(L))
}

# --- Held-out scoring ---------------------------------------------------------------

#' Score a fitted grid on held-out documents (Section 3.7).
#'
#' @param fits named list (by K) of training topicmodels fits.
#' @param foldin_dtm tokens used to infer theta (reconstruction: full eval docs;
#'   completion: fold-in halves).
#' @param score_dtm tokens being scored (reconstruction: same as foldin_dtm;
#'   completion: the complementary halves).
#' @param pi_glob named training baseline (optop_make_baseline(train_dtm)$pi_glob).
#' @param word_at integer vector of K values at which word-level indices are kept.
#' @param min_null null-discrepancy floor delta. NULL = package default
#'   (delta = c); a number sets it independently of the support threshold c --
#'   c governs support coarsening, delta governs which documents enter the
#'   relative-fit estimand, and the two are distinct choices. 0 = strict
#'   positivity. Raw d_model / d_null are returned for every document either
#'   way, so delta can also be varied afterwards without rescoring.
#' @return list(summary, doc, word, minbin_report)
score_heldout <- function(fits, foldin_dtm, score_dtm, pi_glob, c = 1,
                          metrics = c("dev", "chisq", "se"),
                          word_at = NULL, word_metrics = "dev",
                          foldin_seed = 1L, n_threads = 1L, min_null = NULL, resolution = FALSE) {
  K_grid <- as.integer(names(fits))
  theta_list <- lapply(fits, foldin_theta, newdata_dtm = foldin_dtm,
                       seed = foldin_seed)
  phi_list <- lapply(fits, phi_from_fit)

  part <- make_heldout_partition(theta_list, phi_list, score_dtm, pi_glob, c)
  base <- list(pi_glob = pi_glob)
  # held-out baseline word deviance in Poisson form (does not depend on K)
  word_null_dev <- if (!is.null(word_at) && "dev" %in% word_metrics)
    .word_null_dev_poisson(score_dtm, pi_glob) else NULL

  sum_rows <- list(); doc_rows <- list(); word_rows <- list()
  for (k in seq_along(K_grid)) {
    K <- K_grid[k]
    pf <- as_pseudo_fit(theta_list[[k]], phi_list[[k]],
                        doc_ids = rownames(score_dtm),
                        vocab = colnames(score_dtm))
    for (m in metrics) {
      res <- .metric_fun(m)(pf, score_dtm, part, base, macro = TRUE,
                            level = "document", n_threads = n_threads,
                            min_null = min_null)
      sum_rows[[length(sum_rows) + 1L]] <- data.table(
        K = K, metric = m, r2_micro = res$r2, r2_macro = res$r2_macro,
        J_pos = sum(!is.na(res$r2_doc)),
        # docs dropped by the null-discrepancy floor (default min_null = c)
        null_excl_share = res$null_excluded_share %||% NA_real_
      )
      doc_rows[[length(doc_rows) + 1L]] <- .collect_doc(res, K, m, part$L, rownames(score_dtm))
    }
    if (!is.null(word_at) && K %in% word_at) {
      for (m in word_metrics) {
        # Only this guarded branch replaces the invalid package null. Never
        # suppress unrelated warnings or calls whose null is retained.
        resw <- withCallingHandlers(
          .metric_fun(m)(pf, score_dtm, part, base, macro = TRUE,
                         level = "word", n_threads = n_threads),
          optop_word_null_baseline = function(w) {
            if (m == "dev" && !is.null(word_null_dev) &&
                length(word_null_dev) == ncol(score_dtm) &&
                all(is.finite(word_null_dev))) invokeRestart("muffleWarning")
          })
        word_rows[[length(word_rows) + 1L]] <- .collect_word(
          resw, K, m, colnames(score_dtm),
          d_null_override = if (m == "dev") word_null_dev else NULL)
      }
    }
  }

  resolution_rows <- if (resolution) rbindlist(lapply(seq_along(K_grid), function(k)
    support_resolution(theta_list[[k]], phi_list[[k]], score_dtm, part, K_grid[k]))) else NULL
  list(
    resolution = resolution_rows,
    scoring = if (exists("revision_provenance")) revision_provenance() else NULL,
    summary = rbindlist(sum_rows),
    doc = rbindlist(doc_rows),
    word = if (length(word_rows)) rbindlist(word_rows) else NULL,
    minbin_report = part$chisq_min_report,
    partition = part
  )
}

#' Block-wise held-out scoring for very large evaluation sets (E2's conditional
#' truth). The partition rule and every score are per-document, so blocking over
#' documents is exact; Micro/Macro are re-pooled from the stacked doc table.
score_heldout_blocked <- function(fits, foldin_dtm, score_dtm, pi_glob, c = 1,
                                  metrics = "dev", block_docs = 2500L,
                                  foldin_seed = 1L, n_threads = 1L,
                                  min_null = NULL, parallel = FALSE) {
  J <- nrow(score_dtm)
  starts <- seq(1L, J, by = block_docs)
  # blocks are independent (per-document partition and scores); `parallel`
  # spreads them over the current future plan -- only from the MAIN process
  runner <- if (parallel) function(x, f) future.apply::future_lapply(
    x, f, future.seed = NULL) else lapply
  docs <- rbindlist(runner(starts, function(s) {
    rows <- s:min(s + block_docs - 1L, J)
    score_heldout(fits, foldin_dtm[rows, , drop = FALSE],
                  score_dtm[rows, , drop = FALSE], pi_glob, c,
                  metrics = metrics, foldin_seed = foldin_seed,
                  n_threads = n_threads, min_null = min_null)$doc
  }))
  # Pooling uses the SAME retained set as the engine: r2_doc is NA exactly for
  # the documents the floor excluded (delta = min_null, default c; strict
  # positivity at 0), so the floor is never re-derived here from c.
  summary <- docs[, .(
    r2_micro = 1 - sum(d_model[!is.na(r2_doc)]) / sum(d_null[!is.na(r2_doc)]),
    r2_macro = mean(r2_doc, na.rm = TRUE),
    J_pos = sum(!is.na(r2_doc)),
    null_excl_share = mean(is.na(r2_doc))
  ), by = .(K, metric)]
  list(summary = summary, doc = docs)
}

# --- In-sample scoring (native OpTop path) -------------------------------------------

#' @param min_null NULL = package default (null-discrepancy floor at the
#'   partition constant c); a number decouples the floor from the partition
#'   resolution (E5's c-sensitivity arm uses min_null = 1 so the c = 5 support
#'   can be compared at a fixed floor).
score_insample <- function(models, dtm, c = 1, metrics = c("dev", "chisq", "se"),
                           word_at = NULL, word_metrics = "dev",
                           n_threads = 1L, min_null = NULL) {
  K_grid <- as.integer(names(models))
  part <- OpTop::optop_make_partition(models = unname(models), dtm = dtm, c = c,
                                      n_threads = n_threads)
  base <- OpTop::optop_make_baseline(dtm)

  sum_rows <- list(); doc_rows <- list(); word_rows <- list()
  for (k in seq_along(K_grid)) {
    K <- K_grid[k]
    for (m in metrics) {
      res <- .metric_fun(m)(models[[k]], dtm, part, base, macro = TRUE,
                            level = "document", n_threads = n_threads,
                            min_null = min_null)
      sum_rows[[length(sum_rows) + 1L]] <- data.table(
        K = K, metric = m, r2_micro = res$r2, r2_macro = res$r2_macro,
        J_pos = sum(!is.na(res$r2_doc)),
        null_excl_share = res$null_excluded_share %||% NA_real_
      )
      doc_rows[[length(doc_rows) + 1L]] <- .collect_doc(res, K, m, part$L, rownames(dtm))
    }
    if (!is.null(word_at) && K %in% word_at) {
      for (m in word_metrics) {
        resw <- .metric_fun(m)(models[[k]], dtm, part, base, macro = TRUE,
                               level = "word", n_threads = n_threads)
        word_rows[[length(word_rows) + 1L]] <- .collect_word(resw, K, m, colnames(dtm))
      }
    }
  }

  list(
    summary = rbindlist(sum_rows),
    doc = rbindlist(doc_rows),
    word = if (length(word_rows)) rbindlist(word_rows) else NULL,
    minbin_report = part$chisq_min_report,
    partition = part, baseline = base
  )
}

# --- Cross-validation of the local held-out machinery against the package ------------

#' On TRAINING data the local partition builder and the pseudo-fit scoring path
#' must reproduce the native package results exactly (same theta, phi, baseline).
crosscheck_optop <- function(models, dtm, c = 1, tol = 1e-8) {
  theta_list <- lapply(models, theta_from_fit)
  phi_list <- lapply(models, phi_from_fit)
  base <- OpTop::optop_make_baseline(dtm)

  part_pkg <- OpTop::optop_make_partition(unname(models), dtm, c = c)
  part_pkg <- .densify_partition(part_pkg, dtm)   # OpTop >= 0.20 compressed form
  part_loc <- make_heldout_partition(theta_list, phi_list, dtm,
                                     base$pi_glob, c = c)

  checks <- list(
    rare_mask   = identical(unname(part_pkg$rare_mask) > 0,
                            unname(part_loc$rare_mask) > 0),
    L           = isTRUE(all.equal(as.numeric(part_pkg$L),
                                   as.numeric(part_loc$L), tolerance = tol)),
    chisq_min_ok = identical(as.logical(part_pkg$chisq_min_ok),
                             as.logical(part_loc$chisq_min_ok))
  )

  K1 <- names(models)[[length(models)]]
  native <- OpTop::optop_index_deviance(models[[K1]], dtm, part_pkg, base,
                                        macro = TRUE)
  pf <- as_pseudo_fit(theta_list[[K1]], phi_list[[K1]],
                      doc_ids = rownames(dtm), vocab = colnames(dtm))
  pseudo <- OpTop::optop_index_deviance(pf, dtm, part_loc, base, macro = TRUE)

  checks$r2_micro <- isTRUE(all.equal(native$r2, pseudo$r2, tolerance = tol))
  checks$r2_macro <- isTRUE(all.equal(native$r2_macro, pseudo$r2_macro,
                                      tolerance = tol))
  checks$r2_doc <- isTRUE(all.equal(as.numeric(native$r2_doc),
                                    as.numeric(pseudo$r2_doc), tolerance = tol))
  checks
}
