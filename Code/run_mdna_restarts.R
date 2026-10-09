# =============================================================================
# run_mdna_restarts.R -- across-restart dispersion of the MD&A held-out fit.
#
# run_mdna.R fits WarpLDA with n_starts = 3 and keeps the best start; the
# individual restarts are not retained. This script reproduces each production
# restart exactly (the per-start seed in .fit_one_k is
# seed_base + 1000*K + start, so a single-start fit with
# fit_seed_base = 1969 + s reproduces start s of the production run with
# SEED_BASE = 1970), then scores held-out reconstruction and completion at the
# reference topic count, giving the across-restart range of the Deviance Macro
# index that the paper's Section 6.1 sentence reports.
#
# common=2 (default): the three restarts are scored TOGETHER with the production
# candidate grid, on ONE support harmonised over all of them. This is the
# comparison that speaks to the paper's numbers: the restarts are evaluated on
# the cells and retained documents that the reported indices use, enlarged only
# by what the two non-selected restarts add to the harmonisation.
# common=1: support harmonised over the three restarts only (coarser question:
# do the restarts differ on a support of their own?).
# common=0: the pre-revision behaviour -- each restart on its own single-model
# partition, where part of the reported range is a difference of supports.
# Paired per-document differences between restarts are reported with standard
# errors (common >= 1).
#
# Usage: Rscript Code/run_mdna_restarts.R [y1=2015] [y2=2016] [workers=4]
#          [K=50] [sample_n=0] [common=2] [K_grid=10:200:10] [out_suffix=]
# Output: Results/csv/mdna_restart_dispersion_{production_common_|common_|}<TAG>.csv
#         Data/MDNA/mdna_restarts_<TAG>.qs2
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))

.args <- commandArgs(trailingOnly = TRUE)
P <- list(y1 = 2015L, y2 = 2016L, workers = 4L, K = 50L, sample_n = 0L,
          common = 2L, K_grid = "10:200:10", out_suffix = "")
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ",
         paste(names(P), collapse = ", "), call. = FALSE)
  P[[kv[1L]]] <- if (kv[1L] %in% c("out_suffix", "K_grid")) kv[2L]
                 else as.integer(kv[2L])
}
stopifnot(P$common %in% 0:2)
K_PROD <- as.integer(.parse_override_value(P$K_grid))   # production candidate grid
options(optop.revision_suffix = P$out_suffix)
SEED_BASE <- 1970L; C_PART <- 1; N_STARTS_PROD <- 3L
setup_parallel(P$workers)

suffix <- if (P$sample_n > 0L) sprintf("_n%d", P$sample_n) else ""
prep_f <- p_data("MDNA", sprintf("mdna_prep_%d_%d%s.qs2", P$y1, P$y2, suffix))
stopifnot(file.exists(prep_f))
prep <- qs_read(prep_f)
dtm_tr <- prep$dtm_train; dtm_ev <- prep$dtm_ev
TAG <- sprintf("MDNA_%d_%d%s%s", P$y1, P$y2, suffix, P$out_suffix)
log_msg("=== MDNA restarts [%s] K=%d | train %d x %d | held-out %d ===",
        TAG, P$K, nrow(dtm_tr), ncol(dtm_tr), nrow(dtm_ev))

# identical corpus signature to run_mdna.R -> same cache family, but the
# (n_starts = 1, fit_seed_base) spec keys give each restart its own entry
sig <- list(corpus = "mdna_item7_pooled", y1 = P$y1, y2 = P$y2,
            dtm_hash = cfg_hash(list(dim(dtm_tr), Matrix::rowSums(dtm_tr)[1:20],
                                     colnames(dtm_tr)[1:50])))

pi_tr <- OpTop::optop_make_baseline(dtm_tr)$pi_glob
fit_outs <- lapply(seq_len(N_STARTS_PROD), function(s_id) {
  t0 <- proc.time()[["elapsed"]]
  fo <- get_fits_cached(dtm_tr, P$K, "WarpLDA", 1L, 1969L + s_id, sig)
  log_msg("restart %d: fitted/loaded (seed %d) in %.0fs",
          s_id, 1969L + s_id + 1000L * P$K + 1L, proc.time()[["elapsed"]] - t0)
  fo
})
spl <- split_tokens_binomial(dtm_ev, 0.5, seed = SEED_BASE)
score_both <- function(fits) list(
  rec = score_heldout(fits, dtm_ev, dtm_ev, pi_tr, C_PART, "dev",
                      foldin_seed = SEED_BASE),
  com = score_heldout(fits, spl$foldin, spl$score, pi_tr, C_PART, "dev",
                      foldin_seed = SEED_BASE + 1L))
diag_of <- function(s_id, what) {
  v <- fit_outs[[s_id]]$diagnostics[[what]][1L]
  if (is.null(v)) NA_real_ else as.numeric(v)
}

paired <- NULL
if (P$common >= 1L) {
  # pseudo-K labels 1..3: the three fits enter ONE harmonised partition
  fl <- lapply(fit_outs, function(o) o$models[[1L]])
  names(fl) <- as.character(seq_len(N_STARTS_PROD))
  if (P$common == 2L) {
    # production grid (best of three starts per K) joins the harmonisation; the
    # restart labels 1..3 must not collide with a candidate topic count
    stopifnot(min(K_PROD) > N_STARTS_PROD)
    production <- get_fits_cached(dtm_tr, K_PROD, "WarpLDA", 3L, SEED_BASE,
                                  sig)$models
    fl <- c(fl, production)
  }
  t0 <- proc.time()[["elapsed"]]
  sc <- revision_checkpoint(
    "MDNA_restarts",
    list(params = P[setdiff(names(P), "workers")],
         train = digest::digest(dtm_tr, algo = "xxhash64"),
         eval = digest::digest(dtm_ev, algo = "xxhash64")),
    "scores", lapply(score_both(fl), function(x) { x$partition <- NULL; x }))
  log_msg("restarts scored on a common support in %.0fs",
          proc.time()[["elapsed"]] - t0)
  disp <- rbindlist(lapply(names(sc), function(mn) {
    d <- sc[[mn]]$summary[metric == "dev" & K <= N_STARTS_PROD, .(start = K, r2_macro, r2_micro, J_pos)]
    d[, `:=`(K = P$K, eval = mn)]     # := path: "eval" breaks inside .()
    d
  }))
  # is the dispersion larger than evaluation noise? paired, same documents
  paired <- rbindlist(lapply(names(sc), function(mn) {
    g <- paired_gains_all(sc[[mn]]$doc[K <= N_STARTS_PROD], family = "all_pairs")
    g <- g[, .(start_a = K, start_b = K_to, delta_mean, se, n)]
    g[, eval := mn]
    g
  }))
} else {
  # closure arguments must not be named like a column of the table they index
  disp <- rbindlist(lapply(seq_len(N_STARTS_PROD), function(s_id) {
    sc <- score_both(fit_outs[[s_id]]$models)
    rbindlist(lapply(names(sc), function(mn) {
      d <- sc[[mn]]$summary[metric == "dev" & K <= N_STARTS_PROD, .(r2_macro, r2_micro, J_pos, K)]
      d[, `:=`(start = s_id, eval = mn)]
      d
    }))
  }))
}
disp[, `:=`(logLik = vapply(start, diag_of, 0, what = "logLik"),
            perplexity = vapply(start, diag_of, 0, what = "perplexity"))]
setorderv(disp, c("eval", "start"))     # naked NSE breaks on a column named
rng <- disp[, .(min = min(r2_macro),    # "eval"; keep everything quoted
                max = max(r2_macro),
                range = max(r2_macro) - min(r2_macro)), by = "eval"]

f <- proj_path("Results", "csv",
               sprintf("mdna_restart_dispersion_%s%s.csv",
                       if (P$common == 2L) "production_common_" else if (P$common == 1L) "common_" else "", TAG))
fwrite(disp, f)
log_msg("written %s", basename(f))
print(disp); print(rng)
if (!is.null(paired)) {
  fp <- sub("\\.csv$", "_paired.csv", f)
  fwrite(paired, fp); print(paired)
  log_msg("written %s", basename(fp))
}
log_msg(paste0("REPORT: across-restart range of the held-out Deviance Macro ",
               "index at K=%d (%s): reconstruction %.4f, completion %.4f"),
        P$K, if (P$common == 2L) "support common to production grid and restarts"
             else if (P$common == 1L) "support common to the restarts"
             else "single-model supports",
        rng[["range"]][rng[["eval"]] == "rec"],
        rng[["range"]][rng[["eval"]] == "com"])
log_msg("=== MDNA restarts complete ===")

cache_put(list(dispersion = disp, paired = paired, range = rng, params = P,
               support = c("single-model supports", "harmonised over the restarts",
                           "harmonised over the production grid and the restarts")[
                             P$common + 1L]),
          p_data("MDNA", paste0("mdna_restarts_", TAG, ".qs2")), P)
