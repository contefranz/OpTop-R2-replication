# =============================================================================
# run_E6.R -- Word-level dual perspective (Section 3.8.2-3.8.3 + Lemma 2)
#
# The frequency-space twin of E3's document-level Micro-Macro gap. Two
# scenarios share the topic DGP:
#   W1 correctly specified -> at K < K* the retained topics fit common words
#      but miss rare-topic vocabulary, so w-Micro >> w-Macro; the gap shrinks
#      toward zero as K -> K* (underfit-driven divergence, no planting);
#   W2 shared stopword contamination (both sides): a non-topic high-frequency
#      block the model cannot track, leaving a persistent w-Micro - w-Macro gap
#      at K* and placing the planted words in the worst-fit tail (the example
#      Section 3.8.3 names).
#
# Reports, per scenario x seed: word-level w-Micro/w-Macro over the K grid
# (in-sample + held-out reconstruction) under the Section-3.8 filter; the
# per-word R^2 at K* with frequency covariates; and the Lemma-2 residual
# (doc-wise vs word-wise total unbinned fitted deviance, must be ~0).
#
# Usage:  Rscript Code/run_E6.R [smoke|pilot|full] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
rc <- resolve_cfg(cli, "E6", PROFILE)          # honours cfg_from= / out_suffix=
cfg <- rc$cfg; SFX <- rc$out_suffix
if (!rc$from_cache || !is.null(cli$label))
  cfg$label <- run_label(cfg, rc$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E6 [%s] tag=%s S=%d x %d scenarios ===", PROFILE, TAG, cfg$S,
        length(cfg$scenarios))

wf <- cfg$word_filter
units <- CJ(sc_id = names(cfg$scenarios), r = seq_len(cfg$S), sorted = FALSE)

.sim_pair <- function(sc_id, r) {
  sc <- cfg$scenarios[[sc_id]]
  seeds <- make_seeds(cfg$seed_base, r,
                      scenario_id = match(sc_id, names(cfg$scenarios)))
  set.seed(seeds$dgp_seed + 2L)
  Phi <- rdirichlet_mat(cfg$K_true, cfg$W, cfg$beta_DGP)
  contam <- NULL; planted <- integer(0)
  if (isTRUE(sc$contam)) {
    pi_stop <- make_stopword_dist(cfg$W, cfg$n_stop, seed = seeds$dgp_seed + 71L)
    contam <- list(pi_stop = pi_stop, w_mean = sc$w_mean, w_conc = 10,
                   mode = "shared")
    planted <- which(pi_stop > 0)
  }
  sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP,
                           cfg$length_spec, seed = seeds$dgp_seed, Phi = Phi,
                           contamination = contam, doc_prefix = "tr")
  sim_ev <- sim_lda_corpus(cfg$J_eval, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP,
                           cfg$length_spec, seed = seeds$dgp_seed + 7L,
                           Phi = Phi, contamination = contam, doc_prefix = "ev")
  list(seeds = seeds, tr = sim_tr, ev = sim_ev, planted_cols = planted)
}

# ------------------------------ prefit pass -----------------------------------
jobs <- lapply(seq_len(nrow(units)), function(i) {
  sp <- .sim_pair(units$sc_id[i], units$r[i])
  log_seeds("E6", units$sc_id[i], units$r[i], sp$seeds)
  list(dtm = sp$tr$dtm, K_grid = cfg$K_grid, method = cfg$fit_method,
       n_starts = cfg$n_starts, fit_seed_base = sp$seeds$fit_seed_base,
       signature = dgp_signature(cfg, sp$seeds, scenario = units$sc_id[i]))
})
prefit_pool(jobs); rm(jobs); invisible(gc())
log_msg("scoring pass: tail -f %s", proj_path("Data", "scoring.log"))

# ------------------------------ scoring pass ----------------------------------
score_unit <- function(i) {
  revision_checkpoint("E6", cfg, paste0("unit_", i), {
  sc_id <- units$sc_id[i]; r <- units$r[i]
  t0 <- proc.time()[["elapsed"]]
  scoring_log("E6 %s replicate %d/%d: scoring started", sc_id, r, cfg$S)
  sp <- .sim_pair(sc_id, r)
  seeds <- sp$seeds; sim_tr <- sp$tr; sim_ev <- sp$ev

  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base,
                          dgp_signature(cfg, seeds, scenario = sc_id))$models
  base_tr <- OpTop::optop_make_baseline(sim_tr$dtm)

  # word-level indices across the WHOLE grid, in-sample and held-out
  ins <- score_insample(fits, sim_tr$dtm, cfg$c, "dev",
                        word_at = cfg$K_grid, word_metrics = "dev")
  scoring_log("E6 %s r%d: in-sample (word-level) scored (%.0fs)", sc_id, r,
              proc.time()[["elapsed"]] - t0)
  rec <- score_heldout(fits, sim_ev$dtm, sim_ev$dtm, base_tr$pi_glob, cfg$c,
                       "dev", word_at = cfg$K_grid, word_metrics = "dev",
                       foldin_seed = seeds$split_seed)
  scoring_log("E6 %s r%d: held-out (word-level) scored (%.0fs)", sc_id, r,
              proc.time()[["elapsed"]] - t0)

  curve_in <- word_micro_macro_over_k(ins$word, sim_tr$dtm, base_tr$pi_glob,
                                      wf$min_docfreq, wf$min_expected)$curve
  curve_ho <- word_micro_macro_over_k(rec$word, sim_ev$dtm, base_tr$pi_glob,
                                      wf$min_docfreq, wf$min_expected)$curve
  curve <- rbindlist(list(
    copy(curve_in)[, eval := "insample"],
    copy(curve_ho)[, eval := "ho_reconstruction"]
  ))[, `:=`(scenario = sc_id, replicate = r)]

  # per-word R^2 at K* with frequency covariates (held-out), planted flag
  wstar <- word_index_filtered(rec$word[K == cfg$K_true], sim_ev$dtm,
                               base_tr$pi_glob, wf$min_docfreq,
                               wf$min_expected)$word
  wstar[, `:=`(
    log10_docfreq = log10(doc_freq + 1),
    planted = as.integer(sub("word_", "", word_id)) %in% sp$planted_cols,
    scenario = sc_id, replicate = r)]
  word_star <- if (r == 1L)
    wstar[, .(scenario, word_id, r2_word, doc_freq, log10_docfreq, B_w,
              keep, planted)] else NULL

  # Lemma 2 identity on the UNBINNED support (training fit at K*)
  fstar <- fits[[as.character(cfg$K_true)]]
  dw <- ins$word[K == cfg$K_true][match(colnames(sim_tr$dtm), word_id), d_model]
  l2 <- lemma2_residual(theta_from_fit(fstar), phi_from_fit(fstar), sim_tr$dtm,
                         word_d_model = dw)
  stopifnot(abs(l2$rel_resid_pkg) < 1e-10,
            l2$max_word_diff < 1e-8 * max(1, max(abs(dw))))
  lemma2 <- data.table(scenario = sc_id, replicate = r,
                       dev_doc = l2$dev_doc, dev_word = l2$dev_word,
                       resid = l2$resid, rel_resid_pkg = l2$rel_resid_pkg,
                       max_word_diff = l2$max_word_diff, n_floored = l2$n_floored)

  scoring_log("E6 %s replicate %d/%d: DONE in %.0fs", sc_id, r, cfg$S,
              proc.time()[["elapsed"]] - t0)
  log_msg("E6 %s replicate %d/%d scored", sc_id, r, cfg$S)
  list(curve = curve, word_star = word_star, lemma2 = lemma2,
       word_raw = rbindlist(list(copy(ins$word)[, eval := "insample"],
                                 copy(rec$word)[, eval := "ho_reconstruction"]))[
         , `:=`(scenario = sc_id, replicate = r)])
  })
}

res <- future_lapply(seq_len(nrow(units)), score_unit, future.seed = NULL)
pull <- function(name) rbindlist(lapply(res, `[[`, name), fill = TRUE)

out <- list(word_raw = pull("word_raw"), curve = pull("curve"), word_star = pull("word_star"),
            lemma2 = pull("lemma2"), config = cfg)
cache_put(out, p_data("E6", sprintf("e6_results_%s%s.qs2", TAG, SFX)), cfg)
write_result(out$curve, paste0("e6_word_micro_macro", SFX), cfg)
write_result(out$lemma2, paste0("e6_lemma2", SFX), cfg)
log_msg("=== E6 [%s] complete ===", PROFILE)
