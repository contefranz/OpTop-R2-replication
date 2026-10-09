# =============================================================================
# run_E1_selreps.R -- selection distribution on replicated evaluation corpora.
#
# The paper's Table "distribution of K-hat" originally rested on S training
# seeds with ONE evaluation corpus each. This driver keeps the S cached
# training fits (shared fit cache; no refitting) and draws R_sel FRESH
# evaluation corpora per seed (eval seed = dgp_seed + 7 + 1000*rep, disjoint
# from the original +7 corpus), giving S x R_sel selection replicates for the
# eps-adequacy rules and the eval-dependent comparators (held-out perplexity,
# NPMI). optimal_topic is train-only and invariant across evaluation draws,
# so it is not replicated here.
#
# Evaluation-only: fits come from Data/FITS; corpora are regenerated from
# seeds inside workers.
#
# Usage:  Rscript Code/run_E1_selreps.R [smoke|pilot|full] [workers]
#                 [R_sel=10] [key=value ... same overrides as run_E1.R]
# Output: Data/E1/e1_selreps_<TAG>.qs2
#         Results/csv/e1_selreps_khat_<TAG>.csv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
R_SEL <- as.integer(cli$overrides$R_sel %||% "10")
cli$overrides$R_sel <- NULL
stopifnot(R_SEL >= 1L)

rc <- resolve_cfg(cli, "E1", PROFILE)          # honours cfg_from= / out_suffix=
cfg <- rc$cfg; SFX <- rc$out_suffix
if (!rc$from_cache || !is.null(cli$label))
  cfg$label <- run_label(cfg, rc$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E1-selreps [%s] tag=%s S=%d x R_sel=%d ===",
        PROFILE, TAG, cfg$S, R_SEL)

sim_train <- function(seeds) {
  sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                 cfg$length_spec, seed = seeds$dgp_seed, doc_prefix = "tr")
}
sim_eval_rep <- function(seeds, Phi, rep) {
  sim_lda_corpus(cfg$J_eval, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP,
                 cfg$length_spec, seed = seeds$dgp_seed + 7L + 1000L * rep,
                 Phi = Phi, doc_prefix = "ev")
}

# make sure every training fit exists before the scoring fan-out (cache hits
# for the paper configuration; actual fitting only under smoke/pilot configs)
jobs <- list()
for (r in seq_len(cfg$S)) {
  seeds <- make_seeds(cfg$seed_base, r)
  sim_tr <- sim_train(seeds)
  jobs[[length(jobs) + 1L]] <- list(
    dtm = sim_tr$dtm, K_grid = cfg$K_grid, method = cfg$fit_method,
    n_starts = cfg$n_starts, fit_seed_base = seeds$fit_seed_base,
    signature = dgp_signature(cfg, seeds))
  rm(sim_tr)
}
prefit_pool(jobs); rm(jobs); invisible(gc())
log_msg("scoring pass (%d jobs): tail -f %s", cfg$S * R_SEL,
        proj_path("Data", "scoring.log"))

grid <- CJ(r = seq_len(cfg$S), rep = seq_len(R_SEL))

score_job <- function(i) {
  revision_checkpoint("E1_selreps", list(cfg = cfg, R_sel = R_SEL), paste0("unit_", i), {
  r <- grid$r[i]; rep <- grid$rep[i]
  t0 <- proc.time()[["elapsed"]]
  seeds <- make_seeds(cfg$seed_base, r)
  sim_tr <- sim_train(seeds)
  sim_ev <- sim_eval_rep(seeds, sim_tr$Phi, rep)

  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base, dgp_signature(cfg, seeds))$models
  base_tr <- OpTop::optop_make_baseline(sim_tr$dtm)

  fs <- seeds$split_seed + 1000L * rep
  rec <- score_heldout(fits, sim_ev$dtm, sim_ev$dtm, base_tr$pi_glob, cfg$c,
                       "dev", foldin_seed = fs)
  spl <- split_tokens_binomial(sim_ev$dtm, cfg$completion_prop, seed = fs)
  com <- score_heldout(fits, spl$foldin, spl$score, base_tr$pi_glob, cfg$c,
                       "dev", foldin_seed = fs + 1L)

  modes <- list(ho_reconstruction = rec, ho_completion = com)
  khat <- list()
  for (mn in names(modes)) {
    g <- paired_gains(modes[[mn]]$doc, alpha = cfg$sel_alpha)
    for (eps in cfg$eps_grid) {
      khat[[paste(mn, eps)]] <-
        select_k_epsilon(g, eps, cfg$sel_alpha)[
          , `:=`(eval = mn, replicate = r, rep = rep,
                 rule = sprintf("eps_%s", num2tag(eps)))]
    }
  }
  # all-pairs gains + the three rules for every replicate (see run_E1.R), and
  # the per-K summary, whose completion Micro index gives the protocol-matched
  # log-score optimum
  rep_v <- rep      # never write `rep = rep` inside `:=` (column/function clash)
  pairs_all <- rbindlist(lapply(names(modes), function(mn)
    { # no rows when the support collapses (fewer than two retained documents)
      pa <- paired_gains_all(modes[[mn]]$doc, cfg$sel_alpha, NULL, "all_pairs")
      if (nrow(pa)) pa[, `:=`(eval = mn, replicate = r, rep = rep_v)] else NULL
    }), use.names = TRUE)
  sel_all <- rbindlist(lapply(names(modes), function(mn)
    select_k_all_rules(modes[[mn]]$doc, cfg$eps_grid, cfg$sel_alpha)[
      , `:=`(eval = mn, replicate = r, rep = rep_v)]), use.names = TRUE)
  summ <- rbindlist(lapply(names(modes), function(mn)
    copy(modes[[mn]]$summary)[, `:=`(eval = mn, replicate = r, rep = rep_v)]))

  comp <- comparator_metrics(fits, sim_tr$dtm, sim_ev$dtm)
  sel <- select_from_metrics(comp)
  # plain constructor: "eval"/"rep" as literal column names break inside .()
  khat[["comparators"]] <- data.table(
    metric = "comparator", K_hat = sel$K_hat, eps = NA_real_,
    alpha = NA_real_, eval = "comparator", replicate = r, rep = rep,
    rule = sel$rule)

  rep_i <- rep
  qc <- rbindlist(lapply(names(modes), function(mn) {
    d <- modes[[mn]]$summary[metric == "dev" & K == cfg$K_true,
                             .(K, r2_macro)]
    d[, `:=`(eval = mn, replicate = r, rep = rep_i)]
    d
  }))
  scoring_log("E1-selreps r%d rep%d DONE in %.0fs", r, rep,
              proc.time()[["elapsed"]] - t0)
  list(khat = rbindlist(khat, fill = TRUE), qc = qc, pairs_all = pairs_all,
       sel_all = sel_all, summary = summ)
  })
}

res <- future_lapply(seq_len(nrow(grid)), score_job, future.seed = NULL)
khat <- rbindlist(lapply(res, `[[`, "khat"), fill = TRUE)
qc <- rbindlist(lapply(res, `[[`, "qc"), fill = TRUE)

out <- list(khat = khat, qc = qc,
            pairs_all = rbindlist(lapply(res, `[[`, "pairs_all"), use.names = TRUE),
            sel_all = rbindlist(lapply(res, `[[`, "sel_all"), use.names = TRUE),
            summary = rbindlist(lapply(res, `[[`, "summary"), use.names = TRUE),
            R_sel = R_SEL, config = cfg)
cache_put(out, p_data("E1", sprintf("e1_selreps_%s%s.qs2", TAG, SFX)), cfg)
write_result(khat, paste0("e1_selreps_khat", SFX), cfg)
write_result(out$sel_all, paste0("e1_selreps_selection_rules", SFX), cfg)

shares <- khat[rule == "eps_0p01" & metric == "dev",
               .(pr_exact = mean(K_hat == cfg$K_true, na.rm = TRUE),
                 mode_K = as.integer(names(sort(-table(K_hat)))[1L]),
                 mean_K = mean(K_hat, na.rm = TRUE), n = .N),
               by = "eval"]   # quoted: a column literally named "eval"
log_msg("selection shares over %d replicates (dev, eps=0.01):", nrow(grid))
print(shares)
cmp <- khat[eval == "comparator",
            .(pr_exact = mean(K_hat == cfg$K_true, na.rm = TRUE),
              mode_K = as.integer(names(sort(-table(K_hat)))[1L]),
              mean_K = mean(K_hat, na.rm = TRUE), n = .N), by = rule]
print(cmp)
log_msg("=== E1-selreps [%s] complete ===", PROFILE)
