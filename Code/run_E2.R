# =============================================================================
# run_E2.R -- Inference validation via conditional Monte Carlo (Prop 2, Def 1)
#
# The theory is conditional on the training fit, so each training grid is
# fitted ONCE (shared cache with E1) and only evaluation sets are replicated:
#   * conditional truth mu^ho(K), true adjacent gains, true Micro-Macro gap
#     from one very large evaluation set;
#   * R_eval fresh evaluation sets per J_ev: CI coverage for the Macro index,
#     size/power of the paired adjacent-gain test, gap-CI coverage,
#     K_hat selection distribution, t-statistics for QQ plots.
#
# Usage:  Rscript Code/run_E2.R [smoke|pilot|full] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
rc <- resolve_cfg(cli, "E2", PROFILE)          # honours cfg_from= / out_suffix=
cfg <- rc$cfg; SFX <- rc$out_suffix
if (!rc$from_cache || !is.null(cli$label))
  cfg$label <- run_label(cfg, rc$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E2 [%s] tag=%s S_train=%d R=%d ===", PROFILE, TAG,
        cfg$S_train, cfg$R_eval)
log_msg(paste("E2 plan: %d corpora x |K|=%d x %d start(s) = %d LDA fits",
              "(cached; shared with E1). Replication below is eval-only",
              "(fold-in scoring, no refits)."),
        cfg$S_train, length(cfg$K_grid), cfg$n_starts,
        cfg$S_train * length(cfg$K_grid) * cfg$n_starts)
log_msg("per-rep scoring ticks: tail -f %s", proj_path("Data", "scoring.log"))

acc <- list(truth = list(), cover = list(), gains = list(), gap = list(),
            khat = list(), tstats = list())

for (t in seq_len(cfg$S_train)) {
  seeds <- make_seeds(cfg$seed_base, t)          # same schedule as E1 -> cache hits
  log_seeds("E2", "base", t, seeds)

  sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP,
                           cfg$length_spec, seed = seeds$dgp_seed,
                           doc_prefix = "tr")
  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base, dgp_signature(cfg, seeds))$models
  base_tr <- OpTop::optop_make_baseline(sim_tr$dtm)
  pi_tr <- base_tr$pi_glob

  # --- conditional truth from one mega evaluation set -------------------------
  log_msg("E2 corpus %d/%d: conditional truth (fold-in on J=%d eval set)",
          t, cfg$S_train, cfg$J_truth)
  t_tru <- proc.time()[["elapsed"]]
  truth_block <- revision_checkpoint("E2", cfg, paste0("truth_", t), {
  sim_mega <- sim_lda_corpus(cfg$J_truth, cfg$W, cfg$K_true,
                             cfg$alpha_DGP, cfg$beta_DGP, cfg$length_spec,
                             seed = seeds$eval_seed, Phi = sim_tr$Phi,
                             doc_prefix = "mega")
  tru <- score_heldout_blocked(fits, sim_mega$dtm, sim_mega$dtm, pi_tr,
                               cfg$c, metrics = "dev",
                               foldin_seed = seeds$eval_seed, parallel = TRUE)
  log_msg("E2 corpus %d/%d: conditional truth done in %.1f min",
          t, cfg$S_train, (proc.time()[["elapsed"]] - t_tru) / 60)
  # The "truth" is the mean of ONE reference sample of J_truth documents, so it
  # carries sampling error of its own. Its standard errors are stored: coverage
  # of a noisy target is NOT coverage of a known parameter, and the reference
  # precision must be reportable (J_truth = 2000 gives a reference SE of half
  # the evaluation SE at J_ev = 500).
  mu_true <- macro_ci(tru$doc)[metric == "dev", .(K, mu = r2_macro, mu_se = se,
                                                  n_ref = n)]
  gains_true <- paired_gains(tru$doc)[, .(K, delta_true = delta_mean,
                                          delta_true_se = se)]
  gap_true <- micro_macro_gap_ci(tru$doc)[, .(K, gap_true = gap,
                                              gap_true_se = se)]
  list(mu_true = mu_true, gains_true = gains_true, gap_true = gap_true)
  })
  mu_true <- truth_block$mu_true; gains_true <- truth_block$gains_true
  gap_true <- truth_block$gap_true
  acc$truth[[t]] <- cbind(
    train_seed = t,
    Reduce(function(a, b) merge(a, b, by = "K", all = TRUE),
           list(mu_true, gains_true, gap_true))
  )
  invisible(gc())

  # --- replicated evaluation sets ----------------------------------------------
  for (ji in seq_along(cfg$J_ev_grid)) {
    J_ev <- cfg$J_ev_grid[ji]
    log_msg("E2 corpus %d/%d: J_ev=%d, %d eval replications (fold-in scoring, 0 new fits)",
            t, cfg$S_train, J_ev, cfg$R_eval)
    reps <- future_lapply(seq_len(cfg$R_eval), function(rep) {
      revision_checkpoint("E2", cfg, sprintf("t%d_j%d_r%d", t, ji, rep), {
      ev_seed <- seeds$eval_seed + 100000L * ji + 17L * rep
      sim_ev <- sim_lda_corpus(J_ev, cfg$W, cfg$K_true,
                               cfg$alpha_DGP, cfg$beta_DGP,
                               cfg$length_spec, seed = ev_seed,
                               Phi = sim_tr$Phi, doc_prefix = "ev")
      sc <- score_heldout(fits, sim_ev$dtm, sim_ev$dtm, pi_tr, cfg$c,
                          metrics = "dev", foldin_seed = ev_seed)
      ci <- macro_ci(sc$doc)
      g <- paired_gains(sc$doc, alpha = cfg$sel_alpha)
      gp <- micro_macro_gap_ci(sc$doc)
      kh <- rbindlist(lapply(cfg$eps_grid, function(e)
        select_k_epsilon(g, e, cfg$sel_alpha)))
      sa <- select_k_all_rules(sc$doc, cfg$eps_grid, cfg$sel_alpha)
      scoring_log("E2 corpus %d J_ev=%d: rep %d/%d done", t, J_ev, rep,
                  cfg$R_eval)
      list(
        ci = ci[, .(K, r2_macro, se, lwr, upr)][, rep := rep],
        g = g[, .(K, delta_mean, se, ub_onesided)][, rep := rep],
        gp = gp[, .(K, gap, se, lwr, upr)][, rep := rep],
        kh = kh[, rep := rep],
        sa = sa[, rep := rep]
      )
      })
    })
    lab <- function(x) x[, `:=`(train_seed = t, J_ev = J_ev)]
    ci_all <- lab(rbindlist(lapply(reps, `[[`, "ci")))
    g_all  <- lab(rbindlist(lapply(reps, `[[`, "g")))
    gp_all <- lab(rbindlist(lapply(reps, `[[`, "gp")))
    kh_all <- lab(rbindlist(lapply(reps, `[[`, "kh")))
    acc$sel_all[[length(acc$sel_all) + 1L]] <-
      lab(rbindlist(lapply(reps, `[[`, "sa"), use.names = TRUE))

    ci_all <- ci_all[mu_true, on = "K"]
    ci_all[, `:=`(cover = lwr <= mu & mu <= upr, t_stat = (r2_macro - mu) / se)]
    g_all <- g_all[gains_true, on = "K"]
    g_all[, `:=`(
      t_centered = (delta_mean - delta_true) / se,   # size vs conditional truth
      reject_pos = delta_mean / se > qnorm(1 - cfg$sel_alpha)  # power vs 0
    )]
    gp_all <- gp_all[gap_true, on = "K"]
    gp_all[, cover := lwr <= gap_true & gap_true <= upr]

    acc$cover[[length(acc$cover) + 1L]] <- ci_all
    acc$gains[[length(acc$gains) + 1L]] <- g_all
    acc$gap[[length(acc$gap) + 1L]] <- gp_all
    acc$khat[[length(acc$khat) + 1L]] <- kh_all
  }
}

cover <- rbindlist(acc$cover)
gains <- rbindlist(acc$gains)
gap <- rbindlist(acc$gap)
khat <- rbindlist(acc$khat)

cover_tab <- cover[, .(coverage = mean(cover), n = .N), by = .(K, J_ev)]
size_tab <- gains[K >= cfg$K_true,
                  .(size_centered = mean(abs(t_centered) > qnorm(0.975)),
                    n = .N), by = .(K, J_ev)]
power_tab <- gains[K < cfg$K_true,
                   .(power_pos = mean(reject_pos), n = .N), by = .(K, J_ev)]
gap_tab <- gap[, .(coverage = mean(cover), n = .N), by = .(K, J_ev)]
khat_tab <- khat[, .N, by = .(J_ev, eps, K_hat)][order(J_ev, eps, K_hat)]

out <- list(truth = rbindlist(acc$truth), cover = cover, gains = gains,
            gap = gap, khat = khat,
            sel_all = rbindlist(acc$sel_all, use.names = TRUE),
            cover_tab = cover_tab, size_tab = size_tab, power_tab = power_tab,
            gap_tab = gap_tab, khat_tab = khat_tab, config = cfg)
cache_put(out, p_data("E2", sprintf("e2_results_%s%s.qs2", TAG, SFX)), cfg)
write_result(cover_tab, paste0("e2_coverage", SFX), cfg)
write_result(khat_tab, paste0("e2_khat_distribution", SFX), cfg)
log_msg("=== E2 [%s] complete ===", PROFILE)
