# =============================================================================
# stage3_completion_reference.R -- per-fit COMPLETION reference curves (Stage 3 item
# left open on 7 Oct 2026). Mirrors the E2 "conditional truth" computation
# (Code/run_E2.R), but scores the 20,000-document reference corpus of each training
# fit under the completion protocol (binomial half-split, mixture inferred from the
# fold-in half, scoring half scored), with the same floor convention as the E1
# completion indices (package default delta = c = 1). NO MODEL IS FITTED: fits come
# from the shared cache and the run stops if one is missing.
#
# Usage (from the project root):
#   OPTOP_NO_FIT=1 Rscript Code/manuscript/stage3_completion_reference.R smoke   # 1 fit, 300 docs (~1 min)
#   OPTOP_NO_FIT=1 Rscript Code/manuscript/stage3_completion_reference.R pilot   # 1 fit, 20,000 docs (~6 min)
#   OPTOP_NO_FIT=1 Rscript Code/manuscript/stage3_completion_reference.R full    # 10 fits (~1 h on 8 workers)
# Optional second argument: number of workers (default 8).
# Outputs: Results/manuscript/stage3/stage3_completion_reference_curves.csv (per fit x K)
#          Results/manuscript/stage3/stage3_completion_falsecert.csv (E1 completion selections
#          judged against the per-fit curves), Data/E2/e2_completion_truth_rev2[_<mode>].qs2.
# =============================================================================
Sys.setenv(OPTOP_NO_FIT = "1")
suppressPackageStartupMessages(source(here::here("Code", "R", "source_all.R")))
source(here::here("Code", "config", "configs.R"))   # dgp_signature() lives here
options(optop.no_fit = TRUE)
args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "smoke"
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 8L
stopifnot(MODE %in% c("smoke", "pilot", "full"))
SA <- here::here("Results", "manuscript", "stage3"); t0 <- proc.time()[["elapsed"]]
e2 <- qs_read(here::here("Data", "E2", "e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2"))
cfg <- e2$config
n_fits <- if (MODE == "full") cfg$S_train else 1L
J_truth <- if (MODE == "smoke") 300L else cfg$J_truth
sfx <- if (MODE == "full") "" else paste0("_", MODE)
setup_parallel(WORKERS)
log_msg("completion reference: mode=%s fits=%d J_truth=%d K_grid=%s workers=%d", MODE, n_fits, J_truth, paste(range(cfg$K_grid), collapse = "-"), WORKERS)
curves <- list(); docs <- list()
for (t in seq_len(n_fits)) {
  t1 <- proc.time()[["elapsed"]]
  seeds <- make_seeds(cfg$seed_base, t)
  sim_tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP, cfg$length_spec, seed = seeds$dgp_seed, doc_prefix = "tr")
  sig <- dgp_signature(cfg, seeds)
  paths <- vapply(cfg$K_grid, function(K) fit_cache_path_k(sig, K, cfg$fit_method, cfg$n_starts, seeds$fit_seed_base), "")
  stopifnot("a cached fit is missing: refusing to fit" = all(file.exists(paths)))
  fits <- get_fits_cached(sim_tr$dtm, cfg$K_grid, cfg$fit_method, cfg$n_starts, seeds$fit_seed_base, sig)$models
  pi_tr <- OpTop::optop_make_baseline(sim_tr$dtm)$pi_glob
  # the same reference corpus as E2's reconstruction truth (same seed and Phi)
  sim_mega <- sim_lda_corpus(J_truth, cfg$W, cfg$K_true, cfg$alpha_DGP, cfg$beta_DGP, cfg$length_spec, seed = seeds$eval_seed, Phi = sim_tr$Phi, doc_prefix = "mega")
  spl <- split_tokens_binomial(sim_mega$dtm, cfg$completion_prop, seed = seeds$eval_seed + 7L)
  tru <- score_heldout_blocked(fits, spl$foldin, spl$score, pi_tr, cfg$c, metrics = "dev", foldin_seed = seeds$eval_seed + 8L, parallel = TRUE)
  d <- as.data.table(tru$doc)[metric == "dev"]
  mu <- macro_ci(d)[metric == "dev", .(K, mu = r2_macro, mu_se = se, n_ref = n)]
  mu[, `:=`(train_seed = t, J_truth = J_truth, excl_share = d[K == K[1L], mean(is.na(r2_doc))])]
  mu[, A := { m <- mu; sapply(seq_along(m), function(i) if (i < length(m)) max(m[(i + 1):length(m)]) - m[i] else NA_real_) }]
  curves[[t]] <- mu; docs[[t]] <- d[, train_seed := t]
  fwrite(rbindlist(curves), file.path(SA, sprintf("stage3_completion_reference_curves%s.csv", sfx)))
  log_msg("fit %d/%d done in %.1f min: retained %d of %d; mu(40)=%.4f A(40)=%.4f A(50)=%.4f", t, n_fits, (proc.time()[["elapsed"]] - t1) / 60, mu$n_ref[1L], J_truth, mu[K == 40, mu], mu[K == 40, A], mu[K == 50, A])
}
cur <- rbindlist(curves)
qs_save(list(curves = cur, doc = rbindlist(docs), config = cfg, mode = MODE, split_seed_offset = 7L, foldin_seed_offset = 8L), here::here("Data", "E2", sprintf("e2_completion_truth_rev2%s.qs2", sfx)))
cat("\nPer-fit completion remaining gain A_t(K):\n"); print(dcast(cur[K >= 40 & K < 100], K ~ train_seed, value.var = "A"), digits = 3)
if (MODE == "full") {
  sr <- qs_read(here::here("Data", "E1", "e1_selreps_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2"))
  sa <- as.data.table(sr$sel_all)[metric == "dev" & eval == "ho_completion"]
  sa <- merge(sa, cur[, .(replicate = train_seed, K_hat = K, A)], by = c("replicate", "K_hat"), all.x = TRUE)
  res <- sa[, .(n = .N, certified = sum(certified), kstar = sum(certified & K_hat == 40), adequate = sum(certified & A <= eps), false_cert = sum(certified & A > eps), max_A_certified = max(A[certified], na.rm = TRUE)), by = .(rule, eps)][order(rule, eps)]
  cat("\nE1 completion selections (n = 100 per cell) judged against the per-fit completion reference curves:\n"); print(res, digits = 3)
  fwrite(res, file.path(SA, "stage3_completion_falsecert.csv"))
  cat("\nSmallest completion-adequate K per fit:\n"); print(cur[, .(K_0.01 = min(K[!is.na(A) & A <= 0.01]), K_0.005 = min(K[!is.na(A) & A <= 0.005]), A40 = round(A[K == 40], 4)), by = train_seed])
}
log_msg("completion reference complete in %.1f min", (proc.time()[["elapsed"]] - t0) / 60)
