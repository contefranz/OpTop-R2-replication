# =============================================================================
# run_mdna.R -- the real-corpus study: pooled MD&A corpus (fiscal y1 + y2).
#
# Requires Data/MDNA/mdna_prep_<y1>_<y2>[_n<sample>].qs2 (built by prep_mdna.R:
# both fiscal years pooled, one stratified random train/held-out split).
# Reuses the Section-5 machinery wholesale (WarpLDA + shared fit cache, OpTop
# scoring with the null-discrepancy floor, Prop-2 inference, Definition-1
# selection, Section-4 moment tests).
#
# Analyses (paper subsection "An Application to Corporate Disclosures"):
#   A. fit & selection: in-sample / reconstruction / completion over the grid,
#      three families; eps-rule vs perplexity / NPMI / optimal_topic;
#   B. consistency battery at K-hat: Micro/Macro + Prop-2 CI per family, the
#      Micro-Macro gap with its Prop-1(iii) channels, and the null-floor
#      exclusion share (boilerplate index; also by FF12 industry);
#   C. word view: w-Micro/w-Macro over K; worst-fit vocabulary at K-hat;
#   D. Section-4 moment tests on the HELD-OUT set (training instruments,
#      raw statistics + effect sizes; no conditional truth on real data)
#      and the signed residual ranking (over-observed vocabulary).
#
# Fits use ENGINE-DEFAULT priors (the prior-mismatch exercise belongs to the
# simulations, not the application).
#
# Usage:
#   Rscript Code/run_mdna.R [y1=2015] [y2=2016] [workers=4] [K_grid=25:200:25]
#                           [refine_step=10] [refine_span=25] [sample_n=0]
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

# --- CLI ----------------------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
P <- list(y1 = 2015L, y2 = 2016L, workers = 4L, K_grid = "25:200:25",
          refine_step = 10L, refine_span = 25L, sample_n = 0L,
          tests_K = "hat",   # "hat" = K-hat + neighbors; "all" = whole grid
          c_part = 1,        # harmonised-support threshold c (sensitivity runs)
          out_suffix = "")   # appended to the output TAG only; fit cache and
                             # evaluation data are untouched, so sensitivity
                             # runs never overwrite the baseline results
.str_keys <- c("K_grid", "tests_K", "out_suffix")
.num_keys <- c("c_part")
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ",
         paste(names(P), collapse = ", "), call. = FALSE)
  P[[kv[1L]]] <- if (kv[1L] %in% .str_keys) kv[2L]
                 else if (kv[1L] %in% .num_keys) as.numeric(kv[2L])
                 else as.integer(kv[2L])
}
stopifnot(P$tests_K %in% c("hat", "all"))
K_coarse <- .parse_override_value(P$K_grid)
SEED_BASE <- 1970L
EPS_GRID <- c(0.01, 0.005); SEL_ALPHA <- 0.05; C_PART <- P$c_part
stopifnot(is.finite(C_PART), C_PART > 0)
B_STRATA <- 5L; S_STRATA <- 5L; MIN_DOCFREQ <- 5L
setup_parallel(P$workers)

suffix <- if (P$sample_n > 0L) sprintf("_n%d", P$sample_n) else ""
prep_f <- p_data("MDNA", sprintf("mdna_prep_%d_%d%s.qs2", P$y1, P$y2, suffix))
stopifnot(file.exists(prep_f))
prep <- qs_read(prep_f)
dtm_tr <- prep$dtm_train; dtm_ev <- prep$dtm_ev
TAG <- sprintf("MDNA_%d_%d%s%s", P$y1, P$y2, suffix, P$out_suffix)
log_msg("=== MDNA [%s] pooled fiscal %d+%d | train %d x %d | held-out %d ===",
        TAG, P$y1, P$y2, nrow(dtm_tr), ncol(dtm_tr), nrow(dtm_ev))

# corpus signature for the shared fit cache (hash the dtm itself so a rebuilt
# corpus never reuses stale fits)
sig <- list(corpus = "mdna_item7_pooled", y1 = P$y1, y2 = P$y2,
            dtm_hash = cfg_hash(list(dim(dtm_tr), Matrix::rowSums(dtm_tr)[1:20],
                                     colnames(dtm_tr)[1:50])))

fit_grid <- function(K_set) {
  get_fits_cached(dtm_tr, K_set, "WarpLDA", 3L, SEED_BASE, sig)
}

score_all <- function(fits, word_at = NULL) {
  base_tr <- OpTop::optop_make_baseline(dtm_tr)
  pi_tr <- base_tr$pi_glob
  ins <- score_insample(fits, dtm_tr, C_PART, c("dev", "chisq", "se"),
                        word_at = word_at)
  scoring_log("MDNA: in-sample scored (|K|=%d)", length(fits))
  rec <- score_heldout(fits, dtm_ev, dtm_ev, pi_tr, C_PART,
                       c("dev", "chisq", "se"), word_at = word_at,
                       foldin_seed = SEED_BASE)
  scoring_log("MDNA: reconstruction scored")
  spl <- split_tokens_binomial(dtm_ev, 0.5, seed = SEED_BASE)
  com <- score_heldout(fits, spl$foldin, spl$score, pi_tr, C_PART,
                       c("dev", "chisq", "se"), foldin_seed = SEED_BASE + 1L)
  scoring_log("MDNA: completion scored")
  list(ins = ins, rec = rec, com = com, pi_tr = pi_tr)
}

gains_of <- function(sc) {
  rbindlist(lapply(c(rec = "rec", com = "com"), function(mn)
    copy(paired_gains(sc[[mn]]$doc, alpha = SEL_ALPHA))[, eval := mn]),
    fill = TRUE)
}

khat_from_gains <- function(gains_dt) {
  rbindlist(lapply(unique(gains_dt$eval), function(mn) {
    g <- gains_dt[eval == mn]
    rbindlist(lapply(EPS_GRID, function(e)
      select_k_epsilon(g, e, SEL_ALPHA)[
        , `:=`(eval = mn, rule = sprintf("eps_%s", num2tag(e)))]))
  }), fill = TRUE)
}

#' Always-finite K choice: eps-rule (dev, completion, eps = 0.01) -> median of
#' all resolved rules -> gain-elbow heuristic (smallest K after which the mean
#' adjacent dev gain drops below 10% of its maximum; small smoke samples never
#' satisfy the eps bound) -> grid median. Snapped onto K_set.
pick_khat <- function(kh, gains_dt, K_set) {
  k <- kh[eval == "com" & metric == "dev" & rule == "eps_0p01", K_hat][1L]
  if (!is.finite(k))
    k <- suppressWarnings(stats::median(kh$K_hat, na.rm = TRUE))
  if (!is.finite(k)) {
    g <- gains_dt[metric == "dev" & eval == "com"][order(K)]
    if (nrow(g) && any(is.finite(g$delta_mean))) {
      thr <- 0.1 * max(g$delta_mean, na.rm = TRUE)
      idx <- which(g$delta_mean <= thr)
      k <- if (length(idx)) g$K_next[idx[1L]] else max(K_set)
      log_msg("eps-rule unresolved; gain-elbow fallback -> K = %s", k)
    } else {
      k <- stats::median(K_set)
      log_msg("no usable gains; grid-median fallback -> K = %s", k)
    }
  }
  as.integer(K_set[which.min(abs(K_set - k))])
}

# ------------------------- A. fit & selection ----------------------------------
log_msg("A: coarse grid %s", paste(range(K_coarse), collapse = "-"))
fits_c <- fit_grid(K_coarse)
sc_c <- score_all(fits_c$models)
gains_c <- gains_of(sc_c)
kh_c <- khat_from_gains(gains_c)
K_hat_c <- pick_khat(kh_c, gains_c, K_coarse)
log_msg("A: coarse selection K-hat = %s -> refining +/- %d by %d",
        K_hat_c, P$refine_span, P$refine_step)

K_ref <- seq(max(min(K_coarse), K_hat_c - P$refine_span),
             min(max(K_coarse), K_hat_c + P$refine_span), by = P$refine_step)
K_all <- sort(union(K_coarse, K_ref))
fits_out <- fit_grid(K_all)
fits <- fits_out$models

# word-level tracked at every K on the final grid (used by C)
sc <- score_all(fits, word_at = K_all)
gains <- gains_of(sc)
khat <- khat_from_gains(gains)
K_hat <- pick_khat(khat, gains, K_all)
log_msg("A: final selection (dev, completion, eps=0.01): K-hat = %s", K_hat)

comp <- comparator_metrics(fits, dtm_tr, dtm_ev)
sel_comp <- select_from_metrics(comp)
ot <- run_optimal_topic(fits, dtm_tr, alpha = SEL_ALPHA)
sel_comp <- rbindlist(list(sel_comp, data.table(
  rule = if (isTRUE(ot$all_rejected)) "optimal_topic (all rejected; min-stat)"
         else "optimal_topic",
  K_hat = ot$K_hat)))

# ------------------------- B. consistency battery at K-hat ---------------------
# Micro/Macro + Prop-2 CI per family and protocol, at K-hat
battery_ci <- rbindlist(lapply(c(ins = "ins", rec = "rec", com = "com"),
  function(mn) copy(macro_ci(sc[[mn]]$doc[K == K_hat]))[, eval := mn]),
  fill = TRUE)

gap <- rbindlist(lapply(c(ins = "ins", rec = "rec", com = "com"), function(mn)
  copy(micro_macro_gap_ci(sc[[mn]]$doc))[, eval := mn]), fill = TRUE)
decomp <- copy(gap_decomposition(sc$rec$doc[metric == "dev"]))[, eval := "rec"]

# boilerplate index: held-out docs excluded by the null floor at K-hat
# (completion support), overall and by FF12 industry
bp_doc <- sc$com$doc[metric == "dev" & K == K_hat,
                     .(doc_id, d_null, excluded = is.na(r2_doc))]
bp_doc <- merge(bp_doc, prep$dv_ev[, .(doc_id, ff12, sic, fyear)],
                by = "doc_id", all.x = TRUE)
boilerplate <- bp_doc[, .(n = .N, excluded_share = mean(excluded),
                          d_null_median = median(d_null)),
                      by = ff12][order(-excluded_share)]
scoring_log("MDNA: B done (battery CI, gap, decomposition, boilerplate)")

# ------------------------- C. word view ----------------------------------------
wcurve <- word_micro_macro_over_k(sc$rec$word, dtm_ev, sc$pi_tr,
                                  MIN_DOCFREQ, 5)$curve
wstar <- word_index_filtered(sc$rec$word[K == K_hat], dtm_ev, sc$pi_tr,
                             MIN_DOCFREQ, 5)$word
setorder(wstar, r2_word)
scoring_log("MDNA: C done (word curves + worst-fit vocabulary)")

# ------------------------- D. moment tests on the held-out set ------------------
# Section-4 battery with training-built instruments, evaluated on the held-out
# documents at K-hat and its grid neighbors. Raw statistics + effect sizes
# only: on real data there is no conditional truth to center on, so the
# rejection is read jointly with gbar (Remark 8). G accumulated block-wise.
i_hat <- which(K_all == K_hat)
K_test <- if (P$tests_K == "all") {
  K_all   # tests over the whole grid: does more K repair the violations?
} else {
  sort(unique(c(K_hat, K_all[pmax(1L, i_hat - 1L)],
                K_all[pmin(length(K_all), i_hat + 1L)])))
}
log_msg("D: moment tests on held-out set at K in {%s} (tests_K=%s)",
        paste(K_test, collapse = ", "), P$tests_K)

word_tr <- sc$ins$word
Zs_by_K <- lapply(K_test, function(K) {
  ws <- word_tr[K == K & metric == "dev", .(word_id, r2_word)]
  make_instrument_set(dtm_tr, ws, B = B_STRATA, S = S_STRATA,
                      min_docfreq = MIN_DOCFREQ)
})
names(Zs_by_K) <- as.character(K_test)

blocks <- split(seq_len(nrow(dtm_ev)), ceiling(seq_len(nrow(dtm_ev)) / 1500L))
tests_ho <- list(); strata_ho <- list(); resid_sum <- NULL; resid_n <- 0L
for (K in K_test) {
  k <- as.character(K)
  phi_k <- phi_from_fit(fits[[k]])
  Zs <- Zs_by_K[[k]]
  Gs <- lapply(Zs, function(Z) NULL)
  for (bi in seq_along(blocks)) {
    rows <- blocks[[bi]]
    th <- foldin_theta(fits[[k]], dtm_ev[rows, , drop = FALSE],
                       seed = SEED_BASE + bi)
    E <- resid_heldout(th, phi_k, dtm_ev[rows, , drop = FALSE])
    for (nm in names(Zs)) Gs[[nm]] <- rbind(Gs[[nm]], E %*% t(Zs[[nm]]))
    if (K == K_hat) {   # signed residual ranking accumulated at K-hat only
      cs <- colSums(E)
      resid_sum <- if (is.null(resid_sum)) cs else resid_sum + cs
      resid_n <- resid_n + length(rows)
    }
    rm(th, E); invisible(gc(FALSE))
    scoring_log("MDNA: D tests K=%d block %d/%d done", K, bi, length(blocks))
  }
  mt <- lapply(names(Zs), function(nm) moment_test_from_G(Gs[[nm]], nm))
  tests_ho[[k]] <- rbindlist(lapply(mt, `[[`, "result"), fill = TRUE)[, K := K]
  strata_ho[[k]] <- rbindlist(lapply(mt, `[[`, "strata"))[, K := K]
}
tests_ho <- rbindlist(tests_ho)
strata_ho <- rbindlist(strata_ho)

resid_words <- data.table(word = colnames(dtm_ev),
                          resid_mean = resid_sum / resid_n)
freq_share <- function(m) Matrix::colSums(m) / sum(m)
resid_words[, `:=`(share_train = freq_share(dtm_tr)[word],
                   share_ho = freq_share(dtm_ev)[word])]
setorder(resid_words, -resid_mean)
log_msg("D: done; top over-observed word: '%s'", resid_words$word[1L])

# ------------------------- save --------------------------------------------------
out <- list(
  params = P, tag = TAG, K_all = K_all, K_hat = K_hat, K_test = K_test,
  summary = rbindlist(lapply(c(ins = "ins", rec = "rec", com = "com"),
    function(mn) copy(sc[[mn]]$summary)[, eval := mn]), fill = TRUE),
  doc_rec = sc$rec$doc, doc_com = sc$com$doc,
  gains = gains, khat = khat, comparators = comp, sel_comparators = sel_comp,
  optimal_topic = ot[c("K_hat", "all_rejected")],
  battery_ci = battery_ci, gap = gap, decomp = decomp,
  boilerplate = boilerplate, bp_doc = bp_doc,
  wcurve = wcurve, wstar = wstar,
  tests_ho = tests_ho, strata_ho = strata_ho, resid_words = resid_words,
  diagnostics = fits_out$diagnostics, prep_report = prep$report
)
f <- p_data("MDNA", sprintf("mdna_results_%s.qs2", TAG))
cache_put(out, f, list(experiment = "MDNA", profile = "real", tag = TAG))
log_msg("saved %s", basename(f))
log_msg("=== MDNA [%s] complete ===", TAG)
