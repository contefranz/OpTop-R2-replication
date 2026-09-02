# =============================================================================
# run_E4.R -- Moment-based specification tests: size, power, word diagnostics
#
# Size (null DGP): rejection rates of Tests 1-3 at nominal 5%, BOTH raw
#   (H0: mu = 0 -- rejected at scale by any residual imbalance of the fitted
#   model, e.g. VEM smoothing bias; Remark 8) and CENTERED at the conditional
#   truth mu_hat from a large evaluation set (isolates Wald calibration).
#   Effect sizes gbar reported throughout (probability-mass units).
# Power (alternatives; raw H0: mu = 0):
#   * contamination (docvary): doc-specific stopword mixtures -- a genuine
#     misspecification (a SHARED stopword distribution is just one extra LDA
#     topic and gets absorbed; pilot finding);
#   * contamination_eval: clean training, contaminated evaluation documents;
#   * burstiness (Dirichlet-multinomial), drift (eval-side Phi), ctm
#     (logistic-normal theta) -- full profile.
# Companion: R2 across strengths ("fit barely moves, tests fire") and
#   word-level ranks of planted words (eval-only contamination case).
#
# One fold-in per (replication, K): residuals, tests, and R2 share theta.
#
# Usage:  Rscript Code/run_E4.R [smoke|pilot|full] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
cfg <- get_config("E4", PROFILE)
cfg <- apply_overrides(cfg, cli$overrides, strict = FALSE)
cfg$label <- run_label(cfg, cli$overrides, cli$label)
# the battery is evaluated at K_true too, and only at grid members (the null
# branch reuses the full-grid fits, so K_test must be a subset of K_grid)
cfg$K_test <- sort(intersect(unique(c(cfg$K_test, cfg$K_true)), cfg$K_grid))
if (!length(cfg$K_test)) cfg$K_test <- max(cfg$K_grid)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E4 [%s] tag=%s ===", PROFILE, TAG)
log_msg(paste("E4 plan: null arm fits |K|=%d once per corpus (%d corpora);",
              "power arm refits once per contaminated alternative.",
              "Replications are eval-only (fold-in scoring, no refits)."),
        length(cfg$K_grid), cfg$S_train)
log_msg("per-rep scoring ticks: tail -f %s", proj_path("Data", "scoring.log"))

#' Keyed lookup: training word-level deviance scores per tested K (Test 3).
make_word_lookup <- function(word_ins, K_test) {
  out <- data.table(K = K_test)
  out[, ws := lapply(K, function(k)
    word_ins[K == k & metric == "dev", .(word_id, r2_word)])]
  setkey(out, K)
  out
}

#' Instrument sets per tested K (training-only construction).
instruments_by_K <- function(train_dtm, word_tr_l, K_test) {
  out <- lapply(K_test, function(K)
    make_instrument_set(train_dtm, word_tr_l[list(K), ws][[1L]],
                        B = cfg$B_strata, S = cfg$S_strata,
                        min_docfreq = cfg$min_docfreq))
  names(out) <- as.character(K_test)
  out
}

#' Conditional-truth moment centers from one large evaluation set.
estimate_centers <- function(fits_sub, phi_list, Zs_by_K, sim_eval_fun,
                             J_center, seed) {
  ev <- sim_eval_fun(seed, J_center)
  centers <- lapply(names(fits_sub), function(k) {
    th <- foldin_theta(fits_sub[[k]], ev$dtm, seed = seed)
    E <- resid_heldout(th, phi_list[[k]], ev$dtm)
    run_moment_battery(Zs_by_K[[k]], E)$gbars
  })
  names(centers) <- names(fits_sub)
  centers
}

#' Replicated evaluation battery; one fold-in per (rep, K). `label` names the
#' arm in Data/scoring.log (worker-side real-time ticks).
run_battery <- function(fits_sub, phi_list, Zs_by_K, pi_tr, sim_eval_fun,
                        R, K_test, eval_seed0, score_r2 = FALSE,
                        centers = NULL, label = "battery") {
  future_lapply(seq_len(R), function(rep) {
    ev <- sim_eval_fun(eval_seed0 + 17L * rep, cfg$J_ev)
    theta_list <- lapply(fits_sub, foldin_theta, newdata_dtm = ev$dtm,
                         seed = eval_seed0 + rep)
    tests <- rbindlist(lapply(as.character(K_test), function(k) {
      E <- resid_heldout(theta_list[[k]], phi_list[[k]], ev$dtm)
      run_moment_battery(Zs_by_K[[k]], E,
                         centers = centers[[k]])$results[, K := as.integer(k)]
    }), fill = TRUE)
    r2 <- NULL
    if (score_r2) {
      part <- make_heldout_partition(theta_list, phi_list, ev$dtm, pi_tr,
                                     cfg$c)
      r2 <- rbindlist(lapply(as.character(K_test), function(k) {
        pf <- as_pseudo_fit(theta_list[[k]], phi_list[[k]],
                            doc_ids = rownames(ev$dtm),
                            vocab = colnames(ev$dtm))
        res <- OpTop::optop_index_deviance(pf, ev$dtm, part,
                                           list(pi_glob = pi_tr), macro = TRUE)
        data.table(K = as.integer(k), metric = "dev",
                   r2_micro = res$r2, r2_macro = res$r2_macro)
      }))
      r2[, rep := rep]
    }
    scoring_log("E4 %s: rep %d/%d done", label, rep, R)
    list(tests = tests[, rep := rep], r2 = r2)
  })
}

acc <- list(size = list(), power = list(), r2 = list(), strata_ill = NULL,
            word_rank = NULL, word_resid = NULL)

# =============================== SIZE (null DGP) ===============================
for (t in seq_len(cfg$S_train)) {
  seeds <- make_seeds(cfg$seed_base, t)
  log_seeds("E4", "null", t, seeds)
  sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP,
                           cfg$length_spec, seed = seeds$dgp_seed,
                           doc_prefix = "tr")
  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base, dgp_signature(cfg, seeds))$models
  fits_sub <- fits[as.character(cfg$K_test)]
  phi_list <- lapply(fits_sub, phi_from_fit)
  pi_tr <- OpTop::optop_make_baseline(sim_tr$dtm)$pi_glob

  word_ins <- score_insample(fits_sub, sim_tr$dtm, cfg$c, "dev",
                             word_at = cfg$K_test)$word
  word_tr_l <- make_word_lookup(word_ins, cfg$K_test)
  Zs_by_K <- instruments_by_K(sim_tr$dtm, word_tr_l, cfg$K_test)

  sim_eval_null <- function(seed, J) {
    sim_lda_corpus(J, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                   cfg$length_spec, seed = seed, Phi = sim_tr$Phi,
                   doc_prefix = "ev")
  }
  centers <- estimate_centers(fits_sub, phi_list, Zs_by_K, sim_eval_null,
                              cfg$J_center, seeds$eval_seed + 2500000L)
  reps <- run_battery(fits_sub, phi_list, Zs_by_K, pi_tr, sim_eval_null,
                      cfg$R_null, cfg$K_test, seeds$eval_seed + 3000000L,
                      centers = centers, label = sprintf("null t%d", t))
  acc$size[[t]] <- rbindlist(lapply(reps, `[[`, "tests"))[, train_seed := t]
  log_msg("E4 null corpus %d/%d: %d eval replications done (fold-in scoring, 0 new fits)",
          t, cfg$S_train, cfg$R_null)
}

# =============================== POWER (alternatives) ==========================
for (alt in names(cfg$alternatives)) {
  strengths <- cfg$alternatives[[alt]]$strengths
  for (s_val in strengths) {
    for (t in seq_len(cfg$S_train)) {
      # eval-only alternatives keep the CLEAN training corpus (base seeds ->
      # E1/E2 fit-cache reuse); both-sides alternatives get their own seeds.
      eval_only <- alt %in% c("drift", "contamination_eval")
      seeds <- if (eval_only) make_seeds(cfg$seed_base, t)
               else make_seeds(cfg$seed_base, t, scenario_id = 9L)
      pi_stop <- make_stopword_dist(cfg$W, cfg$n_stop,
                                    seed = seeds$dgp_seed + 901L)

      contam_tr <- if (alt == "contamination")
        list(pi_stop = pi_stop, w_mean = s_val, w_conc = 10,
             mode = "docvary", m = 10L) else NULL
      contam_ev <- switch(alt,
        contamination = contam_tr,
        contamination_eval = list(pi_stop = pi_stop, w_mean = s_val,
                                  w_conc = 10, mode = "shared"),
        NULL)
      btau <- if (alt == "burstiness") s_val else Inf
      tsamp <- if (alt == "ctm") theta_logistic_normal(s_val) else NULL
      # document-group x vocabulary-block alternative (both sides): group
      # assignment is regenerated per evaluation set (fresh docs), so only the
      # block layout and weight are shared -- built inside the samplers below.
      gv_tr <- if (alt == "group_vocab")
        make_group_vocab(cfg$J_train, cfg$W, cfg$n_groups, s_val,
                         seeds$dgp_seed + 55L) else NULL

      sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                               cfg$alpha_DGP, cfg$beta_DGP, cfg$length_spec,
                               seed = seeds$dgp_seed,
                               contamination = contam_tr,
                               burst_tau = if (eval_only) Inf else btau,
                               theta_sampler = if (eval_only) NULL else tsamp,
                               group_vocab = gv_tr,
                               doc_prefix = "tr")
      sig <- if (eval_only) dgp_signature(cfg, seeds)
             else dgp_signature(cfg, seeds, scenario = alt,
                                extra = list(strength = s_val))
      K_fit <- if (eval_only) cfg$K_grid else cfg$K_test
      fits <- get_fits_cached(sim_tr$dtm, K_fit, cfg$fit_method, cfg$n_starts,
                              seeds$fit_seed_base, sig)$models
      fits_sub <- fits[as.character(cfg$K_test)]
      phi_list <- lapply(fits_sub, phi_from_fit)
      pi_tr <- OpTop::optop_make_baseline(sim_tr$dtm)$pi_glob

      word_ins <- score_insample(fits_sub, sim_tr$dtm, cfg$c, "dev",
                                 word_at = cfg$K_test)$word
      word_tr_l <- make_word_lookup(word_ins, cfg$K_test)
      Zs_by_K <- instruments_by_K(sim_tr$dtm, word_tr_l, cfg$K_test)

      Phi_ev <- if (alt == "drift")
        drift_phi(sim_tr$Phi, s_val, cfg$beta_DGP, seeds$dgp_seed + 33L)
      else sim_tr$Phi

      sim_eval_alt <- function(seed, J) {
        gv_ev <- if (alt == "group_vocab")
          make_group_vocab(J, cfg$W, cfg$n_groups, s_val, seed + 55L) else NULL
        sim_lda_corpus(J, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                       cfg$length_spec, seed = seed, Phi = Phi_ev,
                       contamination = contam_ev, burst_tau = btau,
                       theta_sampler = tsamp, group_vocab = gv_ev,
                       doc_prefix = "ev")
      }
      reps <- run_battery(fits_sub, phi_list, Zs_by_K, pi_tr, sim_eval_alt,
                          cfg$R_power, cfg$K_test,
                          seeds$eval_seed + 4000000L, score_r2 = TRUE,
                          label = sprintf("power %s=%.3g t%d", alt, s_val, t))
      tests <- rbindlist(lapply(reps, `[[`, "tests"), fill = TRUE)
      r2 <- rbindlist(lapply(reps, `[[`, "r2"))
      acc$power[[length(acc$power) + 1L]] <-
        tests[, `:=`(alt = alt, strength = s_val, train_seed = t)]
      acc$r2[[length(acc$r2) + 1L]] <-
        r2[, `:=`(alt = alt, strength = s_val, train_seed = t)]

      # illustrative rejection + planted-word ranks: eval-only contamination,
      # where the planted words CANNOT be absorbed by the training fit
      if (is.null(acc$strata_ill) && alt == "contamination_eval" &&
          s_val == max(strengths) && t == 1L) {
        ev1 <- sim_eval_alt(seeds$eval_seed + 4000017L, cfg$J_ev)
        th1 <- foldin_theta(fits_sub[[as.character(cfg$K_true)]], ev1$dtm,
                            seed = 1L)
        E1r <- resid_heldout(th1, phi_list[[as.character(cfg$K_true)]],
                             ev1$dtm)
        acc$strata_ill <- run_moment_battery(
          Zs_by_K[[as.character(cfg$K_true)]], E1r)$strata
        scw <- score_heldout(fits_sub, ev1$dtm, ev1$dtm, pi_tr, cfg$c, "dev",
                             word_at = cfg$K_true, foldin_seed = 2L)
        wr <- word_index_filtered(scw$word, ev1$dtm, pi_tr,
                                  cfg$min_docfreq, 5)$word
        planted_ids <- colnames(ev1$dtm)[pi_stop > 0]
        wr[, planted := word_id %in% planted_ids]
        acc$word_rank <- wr
        # Section-4-consistent word diagnostic: mean held-out residual per word
        # (probability-mass units). Systematically over-observed words --
        # planted contamination -- top this ranking mechanically.
        acc$word_resid <- data.table(
          word_id = colnames(ev1$dtm),
          resid_mean = colMeans(E1r),
          planted = colnames(ev1$dtm) %in% planted_ids
        )
      }
      log_msg("E4 power %s=%.3g corpus %d/%d: %d eval replications done (fold-in; training refit once per alternative)",
              alt, s_val, t, cfg$S_train, cfg$R_power)
    }
  }
}

size <- rbindlist(acc$size, fill = TRUE)
power <- rbindlist(acc$power, fill = TRUE)
r2 <- rbindlist(acc$r2)

size_tab <- size[, .(size_raw = mean(pval < 0.05),
                     size_centered = mean(pval_centered < 0.05),
                     gbar_absmax = mean(gbar_absmax), n = .N), by = .(test, K)]
# size-only arms (R_power=0 / alternatives=none) leave power/r2 empty; a
# zero-column rbindlist() result cannot be grouped, so guard before summarising
power_tab <- if (nrow(power)) {
  power[, .(power = mean(pval < 0.05),
            gbar_absmax = mean(gbar_absmax), n = .N),
        by = .(test, K, alt, strength)]
} else power
r2_tab <- if (nrow(r2)) {
  r2[, .(r2_micro = mean(r2_micro), r2_macro = mean(r2_macro)),
     by = .(K, alt, strength)]
} else r2

out <- list(size = size, power = power, r2 = r2,
            size_tab = size_tab, power_tab = power_tab, r2_tab = r2_tab,
            strata_ill = acc$strata_ill, word_rank = acc$word_rank,
            word_resid = acc$word_resid, config = cfg)
cache_put(out, p_data("E4", sprintf("e4_results_%s.qs2", TAG)), cfg)
write_result(size_tab, "e4_size", cfg)
if (nrow(power_tab)) write_result(power_tab, "e4_power", cfg)
log_msg("=== E4 [%s] complete ===", PROFILE)
