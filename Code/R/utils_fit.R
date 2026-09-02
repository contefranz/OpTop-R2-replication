# =============================================================================
# utils_fit.R
# LDA estimation through the NLPstudio::fit_topic_model() facade
# (fit_method: VEM | Gibbs via topicmodels, WarpLDA via text2vec), multi-start
# selection, fold-in of held-out documents, document-completion token splits,
# and construction of the pseudo `nlp_topic_fit` objects that carry
# (theta_hat_eval, Phi_hat_train) through OpTop's official model adapter.
#
# Serialization note: text2vec's WarpLDA R6 object wraps a C++ pointer that
# does NOT survive the qs2 fit cache, so text2vec fits are cached SLIM
# (model_object stripped; dtw/tww/vocab kept) and their held-out fold-in uses
# the engine-agnostic fixed-phi EM below. topicmodels objects serialize fine
# and keep the validated posterior() fold-in. Accessors also accept raw
# topicmodels fits so pre-facade caches remain usable.
# =============================================================================

library(topicmodels)
library(NLPstudio)
library(Matrix)
library(data.table)
library(future.apply)

# --- Engine facade ---------------------------------------------------------------

.engine_spec <- function(fit_method) {
  switch(fit_method,
    VEM     = list(engine = "topicmodels", method = "VEM"),
    Gibbs   = list(engine = "topicmodels", method = "Gibbs"),
    WarpLDA = list(engine = "text2vec",    method = NULL),
    stop("unknown fit_method '", fit_method, "' (use VEM, Gibbs, or WarpLDA)")
  )
}

.is_nlpfit <- function(fit) inherits(fit, "nlp_topic_fit")

#' Strip the non-serializable backend object from text2vec fits before caching.
.slim_fit <- function(fit) {
  if (identical(fit$engine, "text2vec")) fit$model_object <- NULL
  fit
}

# --- Input coercion -------------------------------------------------------------

#' topicmodels-safe input: sparse triplet with dimnames (works for LDA() and
#' posterior() across topicmodels versions; avoids guessing dgCMatrix support).
.as_stm <- function(dtm) {
  if (inherits(dtm, "simple_triplet_matrix")) return(dtm)
  T <- as(as(dtm, "generalMatrix"), "TsparseMatrix")
  slam::simple_triplet_matrix(
    i = T@i + 1L, j = T@j + 1L, v = T@x,
    nrow = nrow(dtm), ncol = ncol(dtm), dimnames = dimnames(dtm)
  )
}

# --- Fitting ---------------------------------------------------------------------

#' Normalize an optional fitted-model prior.  The simulation configuration keeps
#' DGP priors separate from fitting priors: an NA here means "leave the engine at
#' its own default/estimation behaviour" rather than use the DGP value.
.normalize_fit_prior <- function(x, name) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return(NULL)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("%s must be NULL/NA or one finite positive scalar", name),
         call. = FALSE)
  }
  as.numeric(x)
}

#' Fit one K with optional *fitting* priors, distinct from the DGP priors.
#'
#' topicmodels VEM can fix alpha but does not expose a fixed beta/delta value;
#' Gibbs accepts alpha and delta.  WarpLDA accepts both through text2vec's
#' document-topic and topic-word priors.  Unsupported combinations fail rather
#' than silently claiming a prior-mismatch experiment was run.
.fit_one_k <- function(dtm, K, method, n_starts, seed_base,
                       gibbs = list(burnin = 500L, iter = 1000L, thin = 100L),
                       fit_alpha = NA_real_, fit_beta = NA_real_) {
  es <- .engine_spec(method)
  seeds <- as.integer(seed_base + K * 1000L + seq_len(n_starts))
  t0 <- proc.time()[["elapsed"]]
  fit_alpha <- .normalize_fit_prior(fit_alpha, "fit_alpha")
  fit_beta  <- .normalize_fit_prior(fit_beta, "fit_beta")

  if (es$engine == "topicmodels") {
    ctrl_fit <- list(seed = seeds, nstart = n_starts, best = TRUE, verbose = 0L)
    if (es$method == "VEM") {
      if (!is.null(fit_beta)) {
        stop("fit_beta is not supported by topicmodels VEM; use Gibbs or WarpLDA for a fixed topic-word prior",
             call. = FALSE)
      }
      if (!is.null(fit_alpha)) {
        # Without estimate.alpha = FALSE topicmodels treats alpha only as an
        # initializer, which would defeat an explicit prior-mismatch design.
        ctrl_fit$alpha <- fit_alpha
        ctrl_fit$estimate.alpha <- FALSE
      }
    } else {                            # collapsed Gibbs
      ctrl_fit <- c(ctrl_fit, gibbs)
      if (!is.null(fit_alpha)) ctrl_fit$alpha <- fit_alpha
      if (!is.null(fit_beta))  ctrl_fit$delta <- fit_beta
    }
    fit <- NLPstudio::fit_topic_model(dtm, engine = "topicmodels",
                                      model = "lda", k = K, method = es$method,
                                      control = list(fit = ctrl_fit))
    ll <- as.numeric(logLik(fit$model_object))
    px <- NA_real_
    alpha_hat <- mean(fit$model_object@alpha)
  } else {                              # text2vec / WarpLDA: manual multi-start
    fit <- NULL; px <- Inf
    model_control <- list(
      doc_topic_prior = if (is.null(fit_alpha)) 50 / K else fit_alpha,
      topic_word_prior = if (is.null(fit_beta)) 0.1 else fit_beta
    )
    for (s in seeds) {
      set.seed(s)
      f <- NLPstudio::fit_topic_model(
        dtm, engine = "text2vec", model = "lda", k = K,
        control = list(
          model = model_control,
          fit = list(n_iter = 1000L, convergence_tol = 1e-3,
                     progressbar = FALSE)))
      p <- text2vec::perplexity(dtm, topic_word_distribution = f$tww,
                                doc_topic_distribution = f$dtw)
      if (p < px) { px <- p; fit <- f }
    }
    fit <- .slim_fit(fit)
    ll <- NA_real_; alpha_hat <- NA_real_
  }

  list(
    fit = fit,
    diag = data.table(
      K = K, method = method, engine = es$engine, n_starts = n_starts,
      logLik = ll, perplexity = px, alpha_hat = alpha_hat,
      fit_alpha = if (is.null(fit_alpha)) NA_real_ else fit_alpha,
      fit_beta = if (is.null(fit_beta)) NA_real_ else fit_beta,
      elapsed_s = round(proc.time()[["elapsed"]] - t0, 1)
    )
  )
}

#' Fit an LDA grid with multi-start selection (best training log-likelihood
#' for topicmodels engines; lowest training perplexity for WarpLDA).
#'
#' @return list(models = named list by K, diagnostics = data.table)
fit_lda_grid <- function(dtm, K_grid, method = c("VEM", "Gibbs", "WarpLDA"),
                         n_starts = 3L, fit_seed_base = 1L,
                         gibbs = list(burnin = 500L, iter = 1000L, thin = 100L),
                         fit_alpha = NA_real_, fit_beta = NA_real_,
                         parallel = TRUE) {
  method <- match.arg(method)
  runner <- if (parallel) {
    function(ks, f) future_lapply(ks, f, future.seed = NULL)
  } else {
    lapply
  }
  out <- runner(K_grid, function(K)
    .fit_one_k(dtm, K, method, n_starts, fit_seed_base, gibbs,
               fit_alpha = fit_alpha, fit_beta = fit_beta))
  models <- lapply(out, `[[`, "fit")
  names(models) <- as.character(K_grid)
  list(models = models, diagnostics = rbindlist(lapply(out, `[[`, "diag")))
}

# --- Topic-word / doc-topic extraction (engine-agnostic; legacy-compatible) -------

phi_from_fit <- function(fit) {
  ph <- if (.is_nlpfit(fit)) as.matrix(fit$tww)
        else topicmodels::posterior(fit)$terms          # K x W
  ph / rowSums(ph)
}

theta_from_fit <- function(fit) {
  th <- if (.is_nlpfit(fit)) as.matrix(fit$dtw)
        else topicmodels::posterior(fit)$topics         # fitted gamma
  th / rowSums(th)
}

# --- Fold-in ----------------------------------------------------------------------

#' Fixed-phi multinomial-mixture EM fold-in: theta_j <- theta_j * [N_j/(theta_j
#' phi)] phi^T, renormalized. Engine-agnostic; used when the backend cannot
#' fold in after deserialization (cached WarpLDA fits). A handful of sparse
#' matmuls per iteration; per-document EM converges quickly.
foldin_theta_em <- function(phi, foldin_dtm, iter = 200L, tol = 1e-10) {
  J <- nrow(foldin_dtm); K <- nrow(phi)
  Tm <- as(as(foldin_dtm, "generalMatrix"), "TsparseMatrix")
  idx <- cbind(Tm@i + 1L, Tm@j + 1L)
  phiT <- t(phi)
  th <- matrix(1 / K, J, K)
  for (it in seq_len(iter)) {
    P <- th %*% phi
    R <- sparseMatrix(i = Tm@i + 1L, j = Tm@j + 1L,
                      x = Tm@x / pmax(P[idx], 1e-300), dims = dim(foldin_dtm))
    th_new <- th * as.matrix(R %*% phiT)
    th_new <- th_new / rowSums(th_new)
    delta <- max(abs(th_new - th))
    th <- th_new
    if (delta < tol) break
  }
  rownames(th) <- rownames(foldin_dtm)
  th
}

#' Document-topic weights for new documents under training-fitted topics.
#' Dispatch: nlp_topic_fit with a live backend -> predict_topic_model()
#' (posterior E-step for topicmodels, $transform for text2vec); slim fits
#' (cached WarpLDA) -> fixed-phi EM; raw topicmodels fits (legacy caches)
#' -> posterior().
foldin_theta <- function(fit, newdata_dtm, seed = 1L) {
  if (.is_nlpfit(fit)) {
    if (!is.null(fit$model_object)) {
      set.seed(as.integer(seed))
      p <- NLPstudio::predict_topic_model(fit, newdata = newdata_dtm)
      tcols <- grep("^Topic[0-9]+$", names(p), value = TRUE)
      th <- as.matrix(p[, tcols, with = FALSE])
      rownames(th) <- p$doc_id
      th <- th[rownames(newdata_dtm), , drop = FALSE]
    } else {
      th <- foldin_theta_em(phi_from_fit(fit), newdata_dtm)
    }
  } else {
    post <- topicmodels::posterior(fit, newdata = .as_stm(newdata_dtm),
                                   control = list(seed = as.integer(seed),
                                                  verbose = 0L))
    th <- post$topics
  }
  th / rowSums(th)
}

# --- Document-completion token split -----------------------------------------------

#' Binomial 50/50 split of each document's token counts into a fold-in half and
#' a scoring half. Documents with L_j >= 2 are guaranteed a non-empty fold-in
#' side (one token is moved if needed) so that theta can always be inferred;
#' an empty scoring side is legitimate (the document then drops from D_ev,+).
split_tokens_binomial <- function(dtm, prop = 0.5, seed = 1L) {
  set.seed(seed)
  T <- as(as(dtm, "generalMatrix"), "TsparseMatrix")
  x_fold <- rbinom(length(T@x), size = as.integer(T@x), prob = prop)

  dt <- data.table(i = T@i + 1L, j = T@j + 1L, x = as.integer(T@x),
                   xf = as.integer(x_fold))
  dt[, xs := x - xf]

  fix_side <- function(dt, side) {
    # move one token into `side` for docs where that side is empty but L_j >= 2
    tot <- dt[, .(fold = sum(xf), score = sum(xs), L = sum(x)), by = i]
    bad <- tot[L >= 2L & get(side) == 0L, i]
    for (d in bad) {
      other <- if (side == "fold") "xs" else "xf"
      r <- dt[i == d][order(-get(other))][1L]
      sel <- dt[, which(i == d & j == r$j)][1L]
      if (side == "fold") {
        dt[sel, `:=`(xf = xf + 1L, xs = xs - 1L)]
      } else {
        dt[sel, `:=`(xf = xf - 1L, xs = xs + 1L)]
      }
    }
    dt
  }
  dt <- fix_side(dt, "fold")
  dt <- fix_side(dt, "score")

  build <- function(col) {
    keep <- dt[[col]] > 0L
    sparseMatrix(i = dt$i[keep], j = dt$j[keep], x = dt[[col]][keep],
                 dims = dim(dtm), dimnames = dimnames(dtm))
  }
  list(foldin = build("xf"), score = build("xs"))
}

# --- Pseudo model objects for OpTop's adapter ---------------------------------------

#' Wrap (theta, phi) as an `nlp_topic_fit` so OpTop's exported index functions
#' dispatch through their official NLPstudio adapter (utils.R in OpTop):
#' requires numeric matrices with rows summing to 1, unique doc ids and terms.
as_pseudo_fit <- function(theta, phi, doc_ids = rownames(theta),
                          vocab = colnames(phi)) {
  theta <- as.matrix(theta); phi <- as.matrix(phi)
  theta <- theta / rowSums(theta)
  phi   <- phi / rowSums(phi)
  stopifnot(!is.null(doc_ids), !is.null(vocab),
            !anyDuplicated(doc_ids), !anyDuplicated(vocab))
  rownames(theta) <- doc_ids
  colnames(phi)   <- vocab
  structure(
    list(dtw = theta, tww = phi, doc_ids = as.character(doc_ids),
         vocab = as.character(vocab), engine = "optop_sim", model = "pseudo"),
    class = "nlp_topic_fit"
  )
}
