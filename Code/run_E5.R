# =============================================================================
# run_E5.R -- Design sensitivity (appendix): grid extension and threshold c
#
# One corpus (E1 replicate-1 seeds):
#   (a) grid extension (referee MC7): the harmonized support is a union over
#       the K grid, so the SAME model's index changes when the grid grows --
#       score identical fits under nested grids and quantify the drift;
#   (b) c in {1, 5} for Deviance vs Pearson (referee MC8): does c = 5
#       rehabilitate the Pearson index?  Min-bin shares reported throughout.
#
# Usage:  Rscript Code/run_E5.R [smoke|pilot|full] [workers]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

cli <- parse_cli(args)
PROFILE <- cli$profile
cfg <- get_config("E5", PROFILE)
cfg <- apply_overrides(cfg, cli$overrides, strict = FALSE)
cfg$label <- run_label(cfg, cli$overrides, cli$label)
TAG <- run_tag(cfg)
setup_parallel(cli$workers)
log_msg("=== E5 [%s] tag=%s ===", PROFILE, TAG)

seeds <- make_seeds(cfg$seed_base, 1L)
log_seeds("E5", "base", 1L, seeds)

sim <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true,
                      cfg$alpha_DGP, cfg$beta_DGP,
                      cfg$length_spec, seed = seeds$dgp_seed,
                      doc_prefix = "tr")

K_all <- sort(unique(unlist(cfg$grids)))
fits <- get_fits_cached(sim$dtm, K_all, cfg$fit_method, cfg$n_starts,
                        seeds$fit_seed_base,
                        dgp_signature(cfg, seeds, extra = "e5"))$models

# --- (a) grid-extension sensitivity ---------------------------------------------
grid_rows <- list(); minbin_rows <- list()
for (gname in names(cfg$grids)) {
  gK <- cfg$grids[[gname]]
  sc <- score_insample(fits[as.character(gK)], sim$dtm, cfg$c, cfg$metrics)
  grid_rows[[gname]] <- copy(sc$summary)[, grid := gname]
  minbin_rows[[gname]] <- data.table(grid = gname, c = cfg$c,
                                     share = sc$minbin_report$share,
                                     excluded_mass = sc$minbin_report$excluded_mass)
  log_msg("E5 grid '%s' (|K|=%d) scored", gname, length(gK))
}
grid_dt <- rbindlist(grid_rows)

# drift of the reported index for the SAME K as the grid extends
common_K <- Reduce(intersect, cfg$grids)
sens_tab <- dcast(grid_dt[K %in% common_K & metric %in% cfg$metrics],
                  K + metric ~ grid, value.var = "r2_micro")
gnames <- names(cfg$grids)
if (length(gnames) >= 2L) {
  sens_tab[, max_abs_drift :=
    do.call(pmax, c(lapply(gnames[-1], function(g)
      abs(get(g) - get(gnames[1]))), na.rm = TRUE))]
}

# --- (b) c-sensitivity (Deviance vs Pearson) -------------------------------------
main_grid <- cfg$grids[[1L]]
c_rows <- list()
for (cc in cfg$c_grid) {
  # min_null = 1 decouples the null-discrepancy floor from the partition
  # constant: the c-arm compares index LEVELS across partition resolutions at
  # a fixed floor. (Under the coupled default, floor = c = 5 excludes every
  # document at this design point -- reported separately in the paper.)
  sc <- score_insample(fits[as.character(main_grid)], sim$dtm, cc,
                       c("dev", "chisq"), min_null = 1)
  c_rows[[as.character(cc)]] <- copy(sc$summary)[, c_value := cc]
  minbin_rows[[paste0("c", cc)]] <- data.table(
    grid = paste0("main_c", cc), c = cc,
    share = sc$minbin_report$share,
    excluded_mass = sc$minbin_report$excluded_mass)
  log_msg("E5 c=%s sensitivity scored", cc)
}
c_dt <- rbindlist(c_rows)

out <- list(grid_curves = grid_dt, grid_sensitivity = sens_tab,
            c_curves = c_dt, minbin = rbindlist(minbin_rows), config = cfg)
cache_put(out, p_data("E5", sprintf("e5_results_%s.qs2", TAG)), cfg)
write_result(sens_tab, "e5_grid_sensitivity", cfg)
write_result(c_dt, "e5_c_sensitivity", cfg)
log_msg("=== E5 [%s] complete ===", PROFILE)
