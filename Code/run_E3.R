# =============================================================================
# run_E3.R -- Micro-Macro fit heterogeneity (upgraded Experiment 2)
#
# Three scenarios sharing the topic DGP:
#   A heterogeneous lengths (80% Pois(500), 20% Pois(5000));
#   B homogeneous lengths (Pois(1000));
#   C equal lengths, mixed Dirichlet concentration (alpha in {0.2, 2}) --
#     an atypicality-driven gap with NO length heterogeneity (Prop 1iii).
#
# Reported per scenario: held-out + in-sample gap with delta-method CI
# (Remark 6), exact 3-channel decomposition (Prop 1iii), subgroup Macro
# curves, document-level scatters (r_j vs L_j and vs kappa_j), and the
# corrected length summary stats (referee Minor 1).
#
# Usage:  Rscript Code/run_E3.R [smoke|pilot|full] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
cfg <- get_config("E3", PROFILE)
cfg <- apply_overrides(cfg, cli$overrides, strict = FALSE)
cfg$label <- run_label(cfg, cli$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E3 [%s] tag=%s S=%d x %d scenarios ===", PROFILE, TAG, cfg$S,
        length(cfg$scenarios))

.theta_sampler_for <- function(sc) {
  if (identical(sc$theta, "mixed_alpha")) {
    theta_mixed_alpha(sc$alphas, sc$share_first)
  } else NULL
}

units <- CJ(sc_id = names(cfg$scenarios), r = seq_len(cfg$S), sorted = FALSE)

.sim_pair <- function(sc_id, r) {
  sc <- cfg$scenarios[[sc_id]]
  seeds <- make_seeds(cfg$seed_base, r,
                      scenario_id = match(sc_id, names(cfg$scenarios)))
  ts <- .theta_sampler_for(sc)
  sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP, sc$length_spec,
                           seed = seeds$dgp_seed,
                           theta_sampler = ts, doc_prefix = "tr")
  sim_ev <- sim_lda_corpus(cfg$J_eval, cfg$W, cfg$K_true,
                           cfg$alpha_DGP, cfg$beta_DGP, sc$length_spec,
                           seed = seeds$dgp_seed + 7L, Phi = sim_tr$Phi,
                           theta_sampler = ts, doc_prefix = "ev")
  list(seeds = seeds, tr = sim_tr, ev = sim_ev)
}

# ------------------------------ prefit pass -----------------------------------
jobs <- lapply(seq_len(nrow(units)), function(i) {
  sp <- .sim_pair(units$sc_id[i], units$r[i])
  log_seeds("E3", units$sc_id[i], units$r[i], sp$seeds)
  list(dtm = sp$tr$dtm, K_grid = cfg$K_grid, method = cfg$fit_method,
       n_starts = cfg$n_starts, fit_seed_base = sp$seeds$fit_seed_base,
       signature = dgp_signature(cfg, sp$seeds, scenario = units$sc_id[i]))
})
prefit_pool(jobs); rm(jobs); invisible(gc())
log_msg("scoring pass: tail -f %s", proj_path("Data", "scoring.log"))

# ------------------------------ scoring pass ----------------------------------
score_unit <- function(i) {
  sc_id <- units$sc_id[i]; r <- units$r[i]
  t0 <- proc.time()[["elapsed"]]
  scoring_log("E3 %s replicate %d/%d: scoring started", sc_id, r, cfg$S)
  sp <- .sim_pair(sc_id, r)
  seeds <- sp$seeds; sim_tr <- sp$tr; sim_ev <- sp$ev

  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base,
                          dgp_signature(cfg, seeds, scenario = sc_id))$models
  base_tr <- OpTop::optop_make_baseline(sim_tr$dtm)

  ins <- score_insample(fits, sim_tr$dtm, cfg$c, cfg$metrics)
  scoring_log("E3 %s r%d: in-sample scored (%.0fs)", sc_id, r,
              proc.time()[["elapsed"]] - t0)
  rec <- score_heldout(fits, sim_ev$dtm, sim_ev$dtm, base_tr$pi_glob,
                       cfg$c, cfg$metrics, foldin_seed = seeds$split_seed)
  scoring_log("E3 %s r%d: held-out scored (%.0fs)", sc_id, r,
              proc.time()[["elapsed"]] - t0)

  # evaluation-side document groups (length quantiles or theta groups)
  grp_ev <- if (!is.null(sim_ev$doc_group)) {
    data.table(doc_id = rownames(sim_ev$dtm), group = sim_ev$doc_group)
  } else {
    L_ev <- Matrix::rowSums(sim_ev$dtm)
    data.table(doc_id = rownames(sim_ev$dtm),
               group = fifelse(L_ev >= quantile(L_ev, 0.8), "long", "short"))
  }

  modes <- list(insample = ins, ho_reconstruction = rec)
  summary <- rbindlist(lapply(names(modes), function(mn)
    copy(modes[[mn]]$summary)[, `:=`(eval = mn, replicate = r,
                                     scenario = sc_id)]))
  gap <- rbindlist(lapply(names(modes), function(mn)
    micro_macro_gap_ci(modes[[mn]]$doc)[, `:=`(eval = mn, replicate = r,
                                               scenario = sc_id)]))
  decomp <- rbindlist(lapply(names(modes), function(mn)
    gap_decomposition(modes[[mn]]$doc[metric == "dev"])[
      , `:=`(eval = mn, replicate = r, scenario = sc_id)]))

  sub <- merge(rec$doc[metric == "dev"], grp_ev, by = "doc_id")[
    , .(r2_macro_group = mean(r2_doc, na.rm = TRUE), n = .N),
    by = .(K, group)][, `:=`(replicate = r, scenario = sc_id)]

  L_all <- c(Matrix::rowSums(sim_tr$dtm), Matrix::rowSums(sim_ev$dtm))
  lenstats <- data.table(
    scenario = sc_id, replicate = r,
    mean = mean(L_all), sd = sd(L_all),
    p25 = quantile(L_all, .25), median = median(L_all),
    p75 = quantile(L_all, .75)
  )

  scatter <- NULL
  if (r == 1L) {
    d <- merge(rec$doc[metric == "dev" & K == cfg$K_true], grp_ev,
               by = "doc_id")
    d[, kappa := d_null / (2 * L)]
    scatter <- d[, .(scenario = sc_id, doc_id, r2_doc, L, kappa, group)]
  }
  scoring_log("E3 %s replicate %d/%d: DONE in %.0fs", sc_id, r, cfg$S,
              proc.time()[["elapsed"]] - t0)
  log_msg("E3 %s replicate %d/%d scored", sc_id, r, cfg$S)
  list(summary = summary, gap = gap, decomp = decomp, subgroup = sub,
       lenstats = lenstats, scatter = scatter)
}

res <- future_lapply(seq_len(nrow(units)), score_unit, future.seed = NULL)
pull <- function(name) rbindlist(lapply(res, `[[`, name), fill = TRUE)

out <- list(
  summary = pull("summary"), gap = pull("gap"),
  decomp = pull("decomp"), subgroup = pull("subgroup"),
  lenstats = pull("lenstats"), scatter = pull("scatter"),
  config = cfg
)
cache_put(out, p_data("E3", sprintf("e3_results_%s.qs2", TAG)), cfg)
write_result(out$gap, "e3_gap", cfg)
write_result(out$decomp, "e3_gap_decomposition", cfg)
write_result(out$lenstats, "e3_length_stats", cfg)
log_msg("=== E3 [%s] complete ===", PROFILE)
