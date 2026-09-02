# =============================================================================
# run_mdna_restarts.R -- across-restart dispersion of the MD&A held-out fit.
#
# run_mdna.R fits WarpLDA with n_starts = 3 and keeps the best start; the
# individual restarts are not retained. This script reproduces each production
# restart exactly (the per-start seed in .fit_one_k is
# seed_base + 1000*K + start, so a single-start fit with
# fit_seed_base = 1969 + s reproduces start s of the production run with
# SEED_BASE = 1970), then scores held-out reconstruction and completion at the
# selected topic count, giving the across-restart range of the Deviance Macro
# index that the paper's Section 6.1 sentence reports.
#
# Usage: Rscript Code/run_mdna_restarts.R [y1=2015] [y2=2016] [workers=4]
#                                         [K=50] [sample_n=0]
# Output: Results/csv/mdna_restart_dispersion_<TAG>.csv
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

.args <- commandArgs(trailingOnly = TRUE)
P <- list(y1 = 2015L, y2 = 2016L, workers = 4L, K = 50L, sample_n = 0L)
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ",
         paste(names(P), collapse = ", "), call. = FALSE)
  P[[kv[1L]]] <- as.integer(kv[2L])
}
SEED_BASE <- 1970L; C_PART <- 1; N_STARTS_PROD <- 3L
setup_parallel(P$workers)

suffix <- if (P$sample_n > 0L) sprintf("_n%d", P$sample_n) else ""
prep_f <- p_data("MDNA", sprintf("mdna_prep_%d_%d%s.qs2", P$y1, P$y2, suffix))
stopifnot(file.exists(prep_f))
prep <- qs_read(prep_f)
dtm_tr <- prep$dtm_train; dtm_ev <- prep$dtm_ev
TAG <- sprintf("MDNA_%d_%d%s", P$y1, P$y2, suffix)
log_msg("=== MDNA restarts [%s] K=%d | train %d x %d | held-out %d ===",
        TAG, P$K, nrow(dtm_tr), ncol(dtm_tr), nrow(dtm_ev))

# identical corpus signature to run_mdna.R -> same cache family, but the
# (n_starts = 1, fit_seed_base) spec keys give each restart its own entry
sig <- list(corpus = "mdna_item7_pooled", y1 = P$y1, y2 = P$y2,
            dtm_hash = cfg_hash(list(dim(dtm_tr), Matrix::rowSums(dtm_tr)[1:20],
                                     colnames(dtm_tr)[1:50])))

pi_tr <- OpTop::optop_make_baseline(dtm_tr)$pi_glob
rows <- list()
for (s in seq_len(N_STARTS_PROD)) {
  t0 <- proc.time()[["elapsed"]]
  fit_out <- get_fits_cached(dtm_tr, P$K, "WarpLDA", 1L, 1969L + s, sig)
  fits <- fit_out$models
  log_msg("restart %d: fitted/loaded (seed %d) in %.0fs",
          s, 1969L + s + 1000L * P$K + 1L, proc.time()[["elapsed"]] - t0)

  rec <- score_heldout(fits, dtm_ev, dtm_ev, pi_tr, C_PART, "dev",
                       foldin_seed = SEED_BASE)
  spl <- split_tokens_binomial(dtm_ev, 0.5, seed = SEED_BASE)
  com <- score_heldout(fits, spl$foldin, spl$score, pi_tr, C_PART, "dev",
                       foldin_seed = SEED_BASE + 1L)
  for (mn in c(rec = "rec", com = "com")) {
    sc <- if (mn == "rec") rec else com
    d <- sc$summary[metric == "dev" & K == P$K, .(K, r2_macro, r2_micro)]
    d[, `:=`(start = s, eval = mn,   # := path: "eval" breaks inside .()
             logLik = fit_out$diagnostics$logLik[1L],
             perplexity = fit_out$diagnostics$perplexity[1L])]
    rows[[paste(s, mn)]] <- d
  }
  log_msg("restart %d: scored in %.0fs", s, proc.time()[["elapsed"]] - t0)
}
disp <- rbindlist(rows)
setorderv(disp, c("eval", "start"))     # naked NSE breaks on a column named
rng <- disp[, .(min = min(r2_macro),    # "eval"; keep everything quoted
                max = max(r2_macro),
                range = max(r2_macro) - min(r2_macro)), by = "eval"]

f <- proj_path("Results", "csv",
               sprintf("mdna_restart_dispersion_%s.csv", TAG))
fwrite(disp, f)
log_msg("written %s", basename(f))
print(disp); print(rng)
log_msg(paste0("REPORT: across-restart range of the held-out Deviance Macro ",
               "index at K=%d: reconstruction %.4f, completion %.4f"),
        P$K, rng[["range"]][rng[["eval"]] == "rec"],
        rng[["range"]][rng[["eval"]] == "com"])
log_msg("=== MDNA restarts complete ===")
