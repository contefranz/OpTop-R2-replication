# =============================================================================
# run_E1.R -- Held-out fit and topic-number choice (headline experiment)
#
# Per corpus seed: fit the K grid on the training split, score
#   (i) in-sample, (ii) held-out-document reconstruction,
#   (iii) held-out-token completion; compute Prop-2 CIs, paired adjacent
# gains, the Definition-1 selection K_hat_{eps,alpha}; comparators
# (held-out perplexity, NPMI, optimal_topic).
# The former Gibbs / multi-start dispersion robustness arms were removed in the
# 4-study consolidation (MC2/Minor-6 evidence is handled in the response letter).
#
# Structure: one pooled prefit pass over every (corpus x K) job (no per-seed
# barrier; longest K first), then scoring parallelized over replicates.
# Corpora are regenerated from seeds inside workers (deterministic, cheap).
#
# Usage:  Rscript Code/run_E1.R [smoke|pilot|full] [E1|E1b|E1c] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
EXPERIMENT <- cli$experiment %||% "E1"   # E1 | E1b | E1c

cfg <- get_config(EXPERIMENT, PROFILE)
cfg <- apply_overrides(cfg, cli$overrides, strict = FALSE)
cfg$label <- run_label(cfg, cli$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== %s [%s] tag=%s S=%d ===", EXPERIMENT, PROFILE, TAG, cfg$S)

sim_train <- function(seeds) {
  sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                 cfg$length_spec, seed = seeds$dgp_seed, doc_prefix = "tr")
}
sim_eval <- function(seeds, Phi) {
  sim_lda_corpus(cfg$J_eval, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                 cfg$length_spec, seed = seeds$dgp_seed + 7L, Phi = Phi,
                 doc_prefix = "ev")
}

# ------------------------------ prefit pass -----------------------------------
jobs <- list()
for (r in seq_len(cfg$S)) {
  seeds <- make_seeds(cfg$seed_base, r)
  log_seeds(EXPERIMENT, "base", r, seeds)
  sim_tr <- sim_train(seeds)
  jobs[[length(jobs) + 1L]] <- list(
    dtm = sim_tr$dtm, K_grid = cfg$K_grid, method = cfg$fit_method,
    n_starts = cfg$n_starts, fit_seed_base = seeds$fit_seed_base,
    signature = dgp_signature(cfg, seeds))
  rm(sim_tr)
}
prefit_pool(jobs); rm(jobs); invisible(gc())
log_msg("scoring pass: tail -f %s", proj_path("Data", "scoring.log"))

# ------------------------------ scoring pass ----------------------------------
score_one <- function(r) {
  t0 <- proc.time()[["elapsed"]]
  scoring_log("E1 replicate %d/%d: scoring started", r, cfg$S)
  seeds <- make_seeds(cfg$seed_base, r)
  sim_tr <- sim_train(seeds)
  sim_ev <- sim_eval(seeds, sim_tr$Phi)

  fit_out <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                             seeds$fit_seed_base, dgp_signature(cfg, seeds))
  fits <- fit_out$models
  base_tr <- OpTop::optop_make_baseline(sim_tr$dtm)

  ins <- score_insample(fits, sim_tr$dtm, cfg$c, cfg$metrics)
  scoring_log("E1 r%d: in-sample scored (%.0fs)", r,
              proc.time()[["elapsed"]] - t0)
  rec <- score_heldout(fits, sim_ev$dtm, sim_ev$dtm, base_tr$pi_glob, cfg$c,
                       cfg$metrics, foldin_seed = seeds$split_seed)
  scoring_log("E1 r%d: reconstruction scored (%.0fs)", r,
              proc.time()[["elapsed"]] - t0)
  spl <- split_tokens_binomial(sim_ev$dtm, cfg$completion_prop,
                               seed = seeds$split_seed)
  com <- score_heldout(fits, spl$foldin, spl$score, base_tr$pi_glob, cfg$c,
                       cfg$metrics, foldin_seed = seeds$split_seed + 1L)
  scoring_log("E1 r%d: completion scored (%.0fs)", r,
              proc.time()[["elapsed"]] - t0)

  modes <- list(insample = ins, ho_reconstruction = rec, ho_completion = com)
  summary <- rbindlist(lapply(names(modes), function(mn)
    copy(modes[[mn]]$summary)[, `:=`(eval = mn, replicate = r)]))
  ci <- rbindlist(lapply(names(modes), function(mn)
    macro_ci(modes[[mn]]$doc)[, `:=`(eval = mn, replicate = r)]))

  gains <- list(); khat <- list()
  for (mn in c("ho_reconstruction", "ho_completion")) {
    g <- paired_gains(modes[[mn]]$doc, alpha = cfg$sel_alpha)
    gains[[mn]] <- copy(g)[, `:=`(eval = mn, replicate = r)]
    for (eps in cfg$eps_grid) {
      khat[[paste(mn, eps)]] <-
        select_k_epsilon(g, eps, cfg$sel_alpha)[
          , `:=`(eval = mn, replicate = r,
                 rule = sprintf("eps_%s", num2tag(eps)))]
    }
  }

  comp <- comparator_metrics(fits, sim_tr$dtm, sim_ev$dtm)
  sel <- select_from_metrics(comp)
  ot <- run_optimal_topic(fits, sim_tr$dtm, alpha = cfg$sel_alpha)
  ot_rule <- if (isTRUE(ot$all_rejected))
    "optimal_topic (all rejected; min-stat)" else "optimal_topic"
  sel <- rbindlist(list(sel, data.table(rule = ot_rule, K_hat = ot$K_hat)))
  khat[["comparators"]] <- sel[, .(metric = "comparator", K_hat,
                                   eps = NA_real_, alpha = NA_real_,
                                   eval = "comparator", replicate = r, rule)]

  out <- list(
    summary = summary, ci = ci, gains = rbindlist(gains),
    khat = rbindlist(khat, fill = TRUE),
    comp = copy(comp)[, replicate := r],
    diagnostics = fit_out$diagnostics[, replicate := r],
    minbin = data.table(replicate = r,
                        insample_share = ins$minbin_report$share,
                        heldout_share = rec$minbin_report$share)
  )
  if (r == 1L) {
    out$doc_rep1 <- rbindlist(lapply(names(modes), function(mn)
      copy(modes[[mn]]$doc)[, eval := mn]))
  }
  scoring_log("E1 replicate %d/%d: DONE in %.0fs", r, cfg$S,
              proc.time()[["elapsed"]] - t0)
  log_msg("E1 replicate %d/%d scored", r, cfg$S)
  out
}

res <- future_lapply(seq_len(cfg$S), score_one, future.seed = NULL)

pull <- function(name) rbindlist(lapply(res, `[[`, name), fill = TRUE)

out <- list(
  summary = pull("summary"), ci = pull("ci"), gains = pull("gains"),
  khat = pull("khat"), comparators = pull("comp"),
  diagnostics = pull("diagnostics"),
  minbin = pull("minbin"), doc_rep1 = res[[1L]]$doc_rep1, config = cfg
)
cache_put(out, p_data(cfg$experiment,
                      sprintf("%s_results_%s.qs2", tolower(cfg$experiment), TAG)),
          cfg)
write_result(out$summary, sprintf("%s_summary", tolower(cfg$experiment)), cfg)
write_result(out$khat, sprintf("%s_khat", tolower(cfg$experiment)), cfg)
log_msg("=== %s [%s] complete ===", EXPERIMENT, PROFILE)
