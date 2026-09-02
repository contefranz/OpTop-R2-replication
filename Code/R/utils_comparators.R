# =============================================================================
# utils_comparators.R
# Selection-rule comparators for E1 (user-selected set):
#   * held-out perplexity  (NLPstudio::evaluate_topic_model, argmin over K)
#   * NPMI coherence       (NLPstudio::evaluate_topic_model, argmax over K)
#   * optimal_topic()      (OpTop's chi-squared sequential selector, JMLR 2022)
# =============================================================================

library(data.table)

#' Perplexity and NPMI for every fitted K.
#' Fits with a live backend go through NLPstudio::evaluate_topic_model();
#' slim fits (cached WarpLDA, model_object stripped) get held-out perplexity
#' locally from the EM fold-in (same token-NLL convention) and NPMI from tww.
comparator_metrics <- function(models, train_dtm, eval_dtm) {
  rows <- lapply(names(models), function(k) {
    fit <- models[[k]]
    nf <- if (inherits(fit, "nlp_topic_fit")) fit
          else NLPstudio::as_nlp_topic_fit(fit)
    if (!is.null(nf$model_object)) {
      ev <- NLPstudio::evaluate_topic_model(
        nf, training = train_dtm, newdata = eval_dtm,
        metrics = c("coherence_npmi", "held_out_perplexity"),
        level = "aggregate"
      )
      return(data.table(K = as.integer(k), metric = ev$metric,
                        value = ev$value))
    }
    phi <- phi_from_fit(nf)
    th <- foldin_theta_em(phi, eval_dtm)
    P <- th %*% phi
    Tm <- as(as(eval_dtm, "generalMatrix"), "TsparseMatrix")
    nll <- -sum(Tm@x * log(pmax(P[cbind(Tm@i + 1L, Tm@j + 1L)], 1e-12))) /
      sum(Tm@x)
    npmi <- tryCatch(
      NLPstudio::evaluate_topic_model(nf, training = train_dtm,
                                      metrics = "coherence_npmi",
                                      level = "aggregate")$value,
      error = function(e) NA_real_)
    data.table(K = as.integer(k),
               metric = c("held_out_perplexity", "coherence_npmi"),
               value = c(exp(nll), npmi))
  })
  rbindlist(rows)
}

#' K selections implied by the comparator metrics.
select_from_metrics <- function(comp_dt) {
  rbindlist(list(
    comp_dt[metric == "held_out_perplexity",
            .(rule = "perplexity_min", K_hat = K[which.min(value)])],
    comp_dt[metric == "coherence_npmi",
            .(rule = "npmi_max", K_hat = K[which.max(value)])]
  ))
}

#' OpTop's original chi-squared selector on the proportion-weighted training dfm.
#' Sequential logic: smallest K at which adequacy is NOT rejected at `alpha`.
#' When every K is rejected (common at corpus-scale df: the binary-test
#' limitation the paper's Section 1 discusses), fall back to the package's own
#' rule -- the K minimizing the standardized statistic -- and flag it.
run_optimal_topic <- function(models, train_dtm, alpha = 0.05) {
  out <- tryCatch({
    dfm_prop <- quanteda::dfm_weight(quanteda::as.dfm(train_dtm),
                                     scheme = "prop")
    tab <- suppressWarnings(OpTop::optimal_topic(
      topic_models = unname(models), weighted_dfm = dfm_prop,
      alpha = alpha, selection = "sequential",
      do_plot = FALSE, verbose = FALSE
    ))
    tab <- as.data.table(tab)
    accepted <- tab[pval > alpha, topic]
    all_rejected <- length(accepted) == 0L
    khat <- if (all_rejected) tab[which.min(OpTop), topic] else min(accepted)
    list(table = tab, K_hat = as.integer(khat), all_rejected = all_rejected)
  }, error = function(e) {
    warning("optimal_topic failed: ", conditionMessage(e))
    list(table = data.table(), K_hat = NA_integer_, all_rejected = NA)
  })
  out
}
