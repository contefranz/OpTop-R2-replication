# =============================================================================
# postprocess_rev2.R -- manuscript-facing summaries of the revision batch.
#
# Cache-only: reads the result objects at the EXPLICIT paths the batch runner
# recorded (Results/revision<suffix>/inputs.json; or manifest=<json> mapping
# stage keys to paths) and writes one csv per table the manuscript needs
# (Results/csv/rev2_*<suffix>.csv) plus Data/MDNA/revision_postprocess<suffix>.qs2.
# A section whose input is absent is skipped and listed; a section that fails is
# reported and the others still run. No fit is touched, nothing is re-scored.
#
# Sections (stage that produces the input):
#   S1  E1 selection under the three rules, both DGPs; matched optima      (5)
#   S2  E1 10 x 10 selection replicates                                    (8)
#   S3  E2 coverage against the precise reference; selection by J_ev       (4)
#   S4  E4 size / power / planted words (convention established, not guessed) (7)
#   S5  E6 held-out word curves, frequency diagnostic, Lemma S1 check      (6)
#   S6  MD&A word curves, vocabulary, tests (iid + firm cluster), residual
#       mass at every tested K, selections by design (c, grid; delta = 1)  (2, 3)
#   S7  support resolution: E1, MD&A, E5                                   (5, 2, 9)
#   S8  unbinned protocol-matched comparator                               (10)
#   S9  restart dispersion on the production support                       (11)
#
# Usage: Rscript Code/postprocess_rev2.R [suffix=_rev2] [manifest=<json>] [allow_skip=1]
#   A section whose input object is absent makes the run exit with status 1 (so a
#   replication without the deposited result objects cannot report success);
#   allow_skip=1 restores the batch-time behaviour (skip, list, exit 0).
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))

.args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default) {
  v <- sub(paste0("^", key, "="), "", grep(paste0("^", key, "="), .args, value = TRUE))
  if (length(v)) v[[1L]] else default
}
SFX <- arg("suffix", "_rev2")
man_f <- arg("manifest", proj_path("Results", paste0("revision", SFX), "inputs.json"))
if (!file.exists(man_f)) stop("manifest not found: ", man_f, call. = FALSE)
man <- jsonlite::fromJSON(man_f, simplifyVector = FALSE)
paths <- if (!is.null(man$outputs)) man$outputs else man
path_of <- function(key) { p <- unlist(paths[[key]]); if (length(p) == 1L) proj_path(p) else NA_character_ }
EPS <- c(0.01, 0.005); ALPHA <- 0.05
OUT <- list(); skipped <- character(); failed <- character()
emit <- function(dt, name) {
  f <- p_results("csv", sprintf("rev2_%s%s.csv", name, SFX))
  fwrite(dt, f); log_msg("  wrote %s (%d rows)", basename(f), nrow(dt)); invisible(dt)
}
section <- function(id, keys, expr) {
  miss <- keys[!vapply(keys, function(k) isTRUE(file.exists(path_of(k))), NA)]
  if (length(miss)) {
    skipped <<- c(skipped, sprintf("%s (no output of: %s)", id, paste(miss, collapse = ", ")))
    return(invisible())
  }
  log_msg("--- %s", id)
  tryCatch(force(expr), error = function(e) {
    failed <<- c(failed, sprintf("%s: %s", id, conditionMessage(e)))
    log_msg("  FAILED: %s", conditionMessage(e))
  })
}
hit <- function(k, k_true) !is.na(k) & k == k_true     # no selection = a miss
mode_of <- function(v) { v <- v[!is.na(v)]
  if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_ }
dist_of <- function(v) { v <- v[!is.na(v)]
  if (!length(v)) return("")
  paste(sprintf("%d:%d", sort(unique(v)), as.integer(table(v))), collapse = " ") }

sel_summary <- function(sa, K_true, step) {
  sa[, .(n = .N, mode_K = mode_of(K_hat), exact = mean(hit(K_hat, K_true)),
         within_one_step = mean(!is.na(K_hat) & abs(K_hat - K_true) <= step),
         mean_K = mean(K_hat, na.rm = TRUE), none_certified = mean(is.na(K_hat)),
         dist = dist_of(K_hat)), by = c("eval", "metric", "rule", "eps")][
           order(eval, metric, rule, -eps)]
}

# NOTE: a section body is a promise evaluated in the calling (global) frame, so
# results are stored with plain `<-` (a `<<-` there would look past the global env).

# ---- S1. E1: three rules, both generating configurations -------------------------
for (k in c("e1_base", "e1_dgp2")) section(sprintf("S1 %s", k), k, {
  e1 <- qs_read(path_of(k)); cfg <- e1$config; Kt <- cfg$K_true
  step <- min(diff(sort(cfg$K_grid)))
  s <- sel_summary(e1$sel_all, Kt, step)
  emit(s, sprintf("e1_selection_rules_%s", k))
  emit(e1$sel_all[metric == "dev", c("eval", "replicate", "rule", "eps", "K_hat",
                                     "certified", "n_comparisons", "max_ub", "M"),
                  with = FALSE], sprintf("e1_selection_by_seed_%s", k))
  # protocol-matched optimum (completion Micro Deviance == completion log score on
  # the same cells, tokens and retained documents) vs the perplexity comparator
  opt <- e1$summary[metric == "dev", .(K_opt = K[which.max(r2_micro)]),
                    by = c("eval", "replicate")]
  ppx <- e1$comparators[metric == "held_out_perplexity",
                        .(K_min = K[which.min(value)]), by = "replicate"]
  agree <- merge(opt[opt[["eval"]] == "ho_reconstruction", c("replicate", "K_opt"),
                     with = FALSE], ppx, by = "replicate")
  optd <- opt[, .(n = .N, exact = mean(K_opt == Kt), mean_K = mean(K_opt),
                  dist = dist_of(K_opt)), by = "eval"]
  optd[, perplexity_dist := dist_of(ppx$K_min)]
  optd[, seeds_perplexity_eq_micro_rec := sum(agree$K_opt == agree$K_min)]
  emit(optd, sprintf("e1_matched_optimum_%s", k))
  excl <- e1$summary[K == Kt, .(excl_mean = mean(null_excl_share),
                                excl_min = min(null_excl_share),
                                excl_max = max(null_excl_share)), by = c("metric", "eval")]
  emit(excl, sprintf("e1_floor_exclusion_%s", k))
  OUT[[k]] <- list(selection = s, matched = optd, exclusion = excl)
})

# ---- S2. selection on 10 x 10 evaluation corpora -------------------------------------
section("S2 selreps", "selreps", {
  z <- qs_read(path_of("selreps")); cfg <- z$config
  s <- sel_summary(z$sel_all, cfg$K_true, min(diff(sort(cfg$K_grid))))
  emit(s, "selreps_selection_rules"); OUT$selreps <- s
})

# ---- S3. E2: coverage against a precise reference; selection by J_ev ------------------
for (k in c("e2_base", "e2_dgp2")) section(sprintf("S3 %s", k), k, {
  e2 <- qs_read(path_of(k)); cfg <- e2$config; Jt <- cfg$J_truth; Kt <- cfg$K_true
  tr <- e2$truth
  cal <- function(dt, est, tru, tru_se, label) {
    # objects written by the revised driver carry the reference SE in the table
    d <- if (tru_se %in% names(dt)) copy(dt) else
      merge(dt, tr[, c("train_seed", "K", tru_se), with = FALSE], by = c("train_seed", "K"))
    d[, `:=`(cov_raw = lwr <= get(tru) & get(tru) <= upr,
             cov_joint = abs(get(est) - get(tru)) <=
               qnorm(.975) * sqrt(se^2 + get(tru_se)^2),
             ratio = get(tru_se) / se)]
    fit <- d[, .(raw = mean(cov_raw), joint = mean(cov_joint), ratio = mean(ratio)),
             by = c("J_ev", "train_seed")]
    fit[, .(target = label, coverage = mean(raw), mc_se = sd(raw) / sqrt(.N),
            coverage_joint = mean(joint), mc_se_joint = sd(joint) / sqrt(.N),
            min_fit = min(raw), max_fit = max(raw), n_fits = .N,
            ref_se_over_eval_se = mean(ratio),
            expected_if_calibrated = 2 * pnorm(qnorm(.975) / sqrt(1 + J_ev[1L] / Jt)) - 1),
        by = "J_ev"][order(J_ev)]
  }
  cv <- rbindlist(list(cal(e2$cover, "r2_macro", "mu", "mu_se", "average fit (Macro)"),
                       cal(e2$gap, "gap", "gap_true", "gap_true_se", "Micro-Macro gap")))
  cv[, J_truth := Jt]
  emit(cv, sprintf("e2_coverage_%s", k))
  emit(e2$cover_tab, sprintf("e2_coverage_cells_%s", k))
  sj <- e2$sel_all[metric == "dev", .(n = .N, exact = mean(hit(K_hat, Kt)),
                                      mean_K = mean(K_hat, na.rm = TRUE),
                                      none_certified = mean(is.na(K_hat)),
                                      dist = dist_of(K_hat)), by = c("rule", "eps", "J_ev")][
                                        order(rule, -eps, J_ev)]
  emit(sj, sprintf("e2_selection_by_Jev_%s", k))
  OUT[[k]] <- list(coverage = cv, selection_by_Jev = sj)
})

# ---- S4. E4: size, power, planted words ------------------------------------------------
for (k in c("e4_base", "e4_prior")) section(sprintf("S4 %s", k), k, {
  f <- path_of(k); e4 <- qs_read(f)
  emit(e4$size_tab, sprintf("e4_size_%s", k))
  if (!is.null(e4$power_tab) && nrow(e4$power_tab)) emit(e4$power_tab, sprintf("e4_power_%s", k))
  if (!is.null(e4$r2_tab) && nrow(e4$r2_tab)) emit(e4$r2_tab, sprintf("e4_r2_%s", k))
  if (!is.null(e4$word_rank)) {
    conv <- word_null_convention_of(e4, f)       # stamped -> no conversion
    wr <- copy(e4$word_rank)
    if (!identical(conv, WORD_NULL_CONVENTIONS[["poisson"]]))
      stop("planted-word table of a batch output is not in Poisson form: ", conv)
    kk <- wr[keep == TRUE & !is.na(r2_word)][order(r2_word)]
    np <- sum(kk$planted)
    pl <- data.table(convention = conv, n_words = nrow(kk), n_planted = np,
                     precision_at_n = mean(head(kk$planted, np)),
                     recall_top50 = sum(head(kk$planted, 50L)) / np,
                     median_rank_planted = median(which(kk$planted)),
                     min_r2 = min(kk$r2_word))
    emit(pl, sprintf("e4_planted_words_%s", k)); OUT[[k]] <- pl
  }
})

# ---- S5. E6: held-out word curves with the corrected null ---------------------------------
section("S5 e6", "e6", {
  e6 <- qs_read(path_of("e6"))
  cur <- e6$curve[, .(w_micro = mean(r2_micro_word), w_macro = mean(r2_macro_word),
                      gap = mean(gap), gap_sd = sd(gap), n_words = mean(n_words),
                      n_rep = .N), by = c("scenario", "eval", "K")][order(scenario, eval, K)]
  emit(cur, "e6_word_curves")
  emit(e6$lemma2[, .(units = .N, max_abs_rel_resid = max(abs(rel_resid_pkg)),
                     max_word_diff = max(max_word_diff), floored_cells_mean = mean(n_floored)),
                 by = "scenario"], "e6_lemma_crosspath")
  ws <- e6$word_star[keep == TRUE & !is.na(r2_word)]
  ws[, freq_decile := as.integer(cut(rank(doc_freq, ties.method = "first"),
                                     quantile(rank(doc_freq, ties.method = "first"),
                                              0:10 / 10), include.lowest = TRUE,
                                     labels = FALSE)), by = "scenario"]
  fd <- ws[, .(n = .N, doc_freq_median = as.numeric(median(doc_freq)),
               r2_word_mean = mean(r2_word), r2_word_median = median(r2_word),
               share_below_zero = mean(r2_word < 0)), by = c("scenario", "freq_decile")][
                 order(scenario, freq_decile)]
  fd_cor <- ws[, .(spearman_r2_docfreq = cor(r2_word, doc_freq, method = "spearman"),
                   n = .N), by = "scenario"]
  emit(fd, "e6_word_fit_by_frequency_decile"); emit(fd_cor, "e6_word_fit_frequency_correlation")
  OUT$e6 <- list(curves = cur, by_decile = fd, cor = fd_cor)
})

# ---- S6. MD&A -------------------------------------------------------------------------------
wald_cluster <- function(G, cl, label) {
  G <- as.matrix(G); J <- nrow(G); q <- ncol(G); gbar <- colMeans(G)
  S <- rowsum(sweep(G, 2L, gbar), cl); ng <- nrow(S)
  V <- ng / (ng - 1) * crossprod(S) / J^2
  W <- as.numeric(t(gbar) %*% solve(V) %*% gbar)
  data.table(test = label, stat_cluster = W, df = q,
             pval_cluster = pchisq(W, q, lower.tail = FALSE), n_clusters = ng)
}
strata_from_Z <- function(Z) {
  S_ <- nrow(Z) + 1L; s <- rep(NA_integer_, ncol(Z))
  for (b in seq_len(S_ - 1L)) s[Z[b, ] > 0] <- b
  s[Z[1L, ] < 0] <- S_; s
}
section("S6 mdna", "mdna", {
  x <- qs_read(path_of("mdna"))
  prep <- qs_read(p_data("MDNA", sprintf("mdna_prep_%d_%d%s.qs2", x$params$y1, x$params$y2,
                                         if (x$params$sample_n > 0L)
                                           sprintf("_n%d", x$params$sample_n) else "")))
  # the evaluation data must be the data the object was scored on
  stopifnot(identical(digest::digest(prep$dtm_ev, algo = "xxhash64"), x$input_check$dtm_ev_hash))
  emit(x$wcurve[metric == "dev"], "mdna_word_curves")
  # worst-fit vocabulary under the paper's null, at the reference fit and at the
  # total-gain selection
  pi_tr <- OpTop::optop_make_baseline(prep$dtm_train)$pi_glob
  voc <- rbindlist(lapply(unique(c(x$K_hat, intersect(180L, x$K_all))), function(k_v) {
    w <- word_index_filtered(x$word_rec[K == k_v], prep$dtm_ev, pi_tr,
                             x$input_check$min_docfreq, 5)$word
    w <- w[keep == TRUE & !is.na(r2_word)][order(r2_word)]
    top <- head(w[, .(word_id, r2_word, doc_freq, d_model, d_null)], 15L)
    cbind(K_fit = k_v, rank = seq_len(nrow(top)), top,
          n_words = nrow(w), n_below_minus1 = sum(w$r2_word < -1))
  }))
  emit(voc, "mdna_worst_fit_vocabulary")
  # tests: iid and firm-clustered
  cik <- as.character(prep$dv_ev$cik[match(rownames(prep$dtm_ev), prep$dv_ev$doc_id)])
  tc <- rbindlist(lapply(names(x$moments_ho), function(k)
    rbindlist(lapply(names(x$moments_ho[[k]]), function(nm) {
      G <- x$moments_ho[[k]][[nm]]; stopifnot(identical(rownames(G), rownames(prep$dtm_ev)))
      wald_cluster(G, cik, nm)[, K := as.integer(k)]
    }))))
  tests <- merge(x$tests_ho[, .(K, test, stat, df, pval)], tc[, .(K, test, stat_cluster,
                 pval_cluster, n_clusters)], by = c("K", "test"))
  tests[, cluster_over_iid := stat_cluster / stat]
  emit(tests[order(K, test)], "mdna_moment_tests"); emit(x$strata_ho, "mdna_moment_strata")
  # residual mass displaced between groups, at every tested K (mean DOCUMENT
  # probability mass, equal document weights -- not pooled token mass)
  Zs <- make_instruments_by_K(prep$dtm_train, x$word_train, x$K_test,
                              B = x$input_check$B_strata, S = x$input_check$S_strata,
                              min_docfreq = x$input_check$min_docfreq)
  mass <- rbindlist(lapply(as.character(x$K_test), function(k) {
    r <- x$resid_words_by_K[K == as.integer(k)]
    r <- r$resid_mean[match(colnames(prep$dtm_train), r$word)]
    one <- function(tn, lab) {
      st <- strata_from_Z(Zs[[k]][[tn]])
      d <- data.table(stratum = ifelse(is.na(st), 0L, st), r = r)[
        , .(words = .N, mass_pp = 100 * sum(r)), by = "stratum"][order(stratum)]
      d[, `:=`(K = as.integer(k), partition = lab)]
      d[, half_abs_sum_pp := sum(abs(mass_pp[stratum > 0L])) / 2][]
    }
    rbindlist(list(one("T2_freq_strata", "training frequency (Test 2)"),
                   one("T3_fit_strata", "training fit (Test 3)")))
  }))
  emit(mass, "mdna_residual_mass")
  OUT$mdna <- list(wcurve = x$wcurve, vocabulary = voc, tests = tests, mass = mass)
})
section("S6b mdna designs", c("mdna", "mdna_c05", "mdna_c2", "mdna_g100"), {
  labs <- c(mdna = "c = 1 (baseline)", mdna_c05 = "c = 0.5", mdna_c2 = "c = 2",
            mdna_g100 = "c = 1, grid 10-100")
  des <- rbindlist(lapply(names(labs), function(k) {
    y <- qs_read(path_of(k))
    s <- y$sel_all[metric == "dev", c("eval", "rule", "eps", "K_hat", "certified", "max_ub", "M"),
                   with = FALSE]
    fit50 <- y$summary[metric == "dev" & K == 50L, c("eval", "r2_micro", "r2_macro",
                                                     "null_excl_share"), with = FALSE]
    s <- merge(s, fit50, by = "eval", all.x = TRUE)
    s[, `:=`(design = labs[[k]], c = y$input_check$c_part, delta = y$input_check$min_null,
             n_models = length(y$K_all))][]
  }))
  setcolorder(des, c("design", "c", "delta", "n_models"))
  emit(des[order(design, eval, -eps, rule)], "mdna_design_sensitivity"); OUT$mdna_designs <- des
})

# ---- S7. support resolution ---------------------------------------------------------------------
res_floor <- function(summary_dt, by_) summary_dt[metric == "dev", c(by_, "null_excl_share"), with = FALSE]
for (k in c("e1_base", "e1_dgp2")) section(sprintf("S7 resolution %s", k), k, {
  e1 <- qs_read(path_of(k))
  r <- summarise_resolution(e1$resolution, by = c("eval", "replicate", "K"))
  r <- merge(r, res_floor(e1$summary, c("eval", "replicate", "K")), by = c("eval", "replicate", "K"))
  num <- setdiff(names(r), c("eval", "replicate", "K"))
  emit(r[, lapply(.SD, mean), by = c("eval", "K"), .SDcols = num][order(eval, K)],
       sprintf("support_resolution_%s", k))
})
section("S7 resolution mdna", "mdna", {
  x <- qs_read(path_of("mdna"))
  r <- summarise_resolution(x$resolution, by = c("eval", "K"))
  fl <- x$summary[metric == "dev" & eval %in% c("rec", "com"), c("eval", "K", "null_excl_share"),
                  with = FALSE]
  emit(merge(r, fl, by = c("eval", "K"))[order(eval, K)], "support_resolution_mdna")
})
section("S7 resolution e5", "e5_resolution", {
  z <- qs_read(path_of("e5_resolution"))
  fl <- z$floor[metric == "dev", c("design", "K", "null_excl_share"), with = FALSE]
  emit(merge(z$summary, fl, by = c("design", "K"))[order(design, K)], "support_resolution_e5")
  emit(z$minbin, "support_resolution_e5_pearson_minbin")
})

# ---- S8. unbinned comparator ------------------------------------------------------------------------
section("S8 unbinned", "unbinned", {
  u <- qs_read(path_of("unbinned")); s <- u$summary
  emit(s[order(source, protocol, floor, K)], "unbinned_comparator")
  best <- s[, .(K_best_log_score = K[which.max(log_score_per_token)],
                K_best_micro = K[which.max(r2_micro)], K_best_macro = K[which.max(r2_macro)],
                max_token_share_floor = max(token_share_floor)),
            by = c("source", "protocol", "floor")][order(source, protocol, floor)]
  spread <- s[, .(log_score_range_across_floors = max(log_score_per_token) -
                    min(log_score_per_token)), by = c("source", "protocol", "K")][
                      , .(max_range = max(log_score_range_across_floors)),
                      by = c("source", "protocol")]
  emit(best, "unbinned_comparator_optima"); emit(spread, "unbinned_comparator_floor_sensitivity")
  OUT$unbinned <- list(optima = best, floor_sensitivity = spread)
})

# ---- S9. restarts ---------------------------------------------------------------------------------------
section("S9 restarts", "restarts", {
  z <- qs_read(path_of("restarts"))
  emit(z$dispersion, "mdna_restart_dispersion"); emit(z$range, "mdna_restart_range")
  if (!is.null(z$paired)) emit(z$paired, "mdna_restart_paired")
  OUT$restarts <- z[c("dispersion", "range", "paired", "support")]
})

f_out <- p_data("MDNA", sprintf("revision_postprocess%s.qs2", SFX))
cache_put(OUT, f_out, list(experiment = "revision", suffix = SFX))
log_msg("saved %s [%s]", basename(f_out), paste(names(OUT), collapse = ", "))
if (length(skipped)) log_msg("SKIPPED (input not produced yet): %s", paste(skipped, collapse = " | "))
if (length(skipped) && !identical(arg("allow_skip", "0"), "1")) {
  log_msg("missing inputs: fetch the result objects (./reproduce.sh fetch results) or pass allow_skip=1")
  quit(status = 1L)
}
if (length(failed)) { log_msg("FAILED: %s", paste(failed, collapse = " | ")); quit(status = 1L) }
log_msg("=== postprocess_rev2 complete ===")
