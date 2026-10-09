# Stage 3: MD&A moment tests on COMPLETION residuals at K = 50 and 180 (no refitting; Table S25).
# Needs the deposited fit cache (Data/FITS). ~4 min on 8 workers. Outputs: Results/manuscript/stage3/
# stage3_mdna_completion_{tests,strata,masses}.csv and Data/MDNA/stage3_mdna_completion_G.qs2.
Sys.setenv(OPTOP_NO_FIT = "1")
suppressPackageStartupMessages(source(here::here("Code", "R", "source_all.R")))
options(optop.no_fit = TRUE)
SP <- here::here("Results", "manuscript", "stage3"); t0 <- proc.time()[["elapsed"]]
prep <- qs_read(p_data("MDNA", "mdna_prep_2015_2016.qs2")); dtm_tr <- prep$dtm_train; dtm_ev <- prep$dtm_ev
x <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016_rev2.qs2"))
stopifnot(identical(digest::digest(dtm_ev, algo = "xxhash64"), x$input_check$dtm_ev_hash))
SEED_BASE <- 1970L; K_test <- c(50L, 180L)
sig <- list(corpus = "mdna_item7_pooled", y1 = 2015L, y2 = 2016L,
            dtm_hash = cfg_hash(list(dim(dtm_tr), Matrix::rowSums(dtm_tr)[1:20], colnames(dtm_tr)[1:50])))
setup_parallel(8L)
fits <- get_fits_cached(dtm_tr, K_test, "WarpLDA", 3L, SEED_BASE, sig)$models
cat(sprintf("[%.0fs] fits loaded from cache: %s\n", proc.time()[["elapsed"]] - t0, paste(names(fits), collapse = ",")))
word_tr <- as.data.table(x$word_train)
Zs_by_K <- make_instruments_by_K(dtm_tr, word_tr, K_test, B = 5L, S = 5L, min_docfreq = 5L)
cik <- as.character(prep$dv_ev$cik[match(rownames(dtm_ev), prep$dv_ev$doc_id)])
wald_cluster <- function(G, cl) { G <- as.matrix(G); J <- nrow(G); q <- ncol(G); gbar <- colMeans(G)
  S <- rowsum(sweep(G, 2L, gbar), cl); ng <- nrow(S); V <- ng / (ng - 1) * crossprod(S) / J^2
  W <- as.numeric(t(gbar) %*% solve(V) %*% gbar); list(W = W, p = pchisq(W, q, lower.tail = FALSE), ng = ng, t_cl = gbar / sqrt(diag(V))) }
strata_from_Z <- function(Z) { S_ <- nrow(Z) + 1L; s <- rep(NA_integer_, ncol(Z)); for (b in seq_len(S_ - 1L)) s[Z[b, ] > 0] <- b; s[Z[1L, ] < 0] <- S_; s }
moments <- function(k, foldin_dtm, score_dtm, seed_off, label) {
  keep <- which(Matrix::rowSums(score_dtm) > 0)
  blocks <- split(keep, ceiling(seq_along(keep) / 1500L))
  phi_k <- phi_from_fit(fits[[as.character(k)]]); Zs <- Zs_by_K[[as.character(k)]]
  Gs <- lapply(Zs, function(Z) NULL); mass <- lapply(Zs, function(Z) NULL)
  for (bi in seq_along(blocks)) { rows <- blocks[[bi]]
    th <- foldin_theta(fits[[as.character(k)]], foldin_dtm[rows, , drop = FALSE], seed = SEED_BASE + seed_off + bi)
    E <- resid_heldout(th, phi_k, score_dtm[rows, , drop = FALSE])
    for (nm in names(Zs)) { Gs[[nm]] <- rbind(Gs[[nm]], E %*% t(Zs[[nm]]))
      s <- strata_from_Z(Zs[[nm]]); M <- sapply(sort(unique(na.omit(s))), function(b) rowSums(E[, which(s == b), drop = FALSE]))
      mass[[nm]] <- rbind(mass[[nm]], M) }
    rm(th, E); invisible(gc(FALSE)); cat(sprintf("[%.0fs] %s K=%d block %d/%d (%d docs)\n", proc.time()[["elapsed"]] - t0, label, k, bi, length(blocks), length(rows))) }
  cl <- cik[keep]
  res <- rbindlist(lapply(names(Gs), function(nm) { G <- Gs[[nm]]; rownames(G) <- rownames(dtm_ev)[keep]
    mt <- moment_test_from_G(G, nm); wc <- wald_cluster(G, cl)
    data.table(protocol = label, K = k, test = nm, J = nrow(G), W_iid = mt$result$stat, p_iid = mt$result$pval, df = mt$result$df, W_cluster = wc$W, p_cluster = wc$p, n_clusters = wc$ng) }))
  strata <- rbindlist(lapply(names(Gs), function(nm) { G <- Gs[[nm]]; wc <- wald_cluster(G, cl); mt <- moment_test_from_G(G, nm)
    cbind(protocol = label, K = k, mt$strata, t_cluster = wc$t_cl) }))
  masses <- rbindlist(lapply(names(mass), function(nm) { M <- mass[[nm]]; mb <- colMeans(M) * 100
    Sg <- rowsum(sweep(M, 2L, colMeans(M)), cl); V <- nrow(Sg) / (nrow(Sg) - 1) * crossprod(Sg) / nrow(M)^2
    data.table(protocol = label, K = k, test = nm, group = seq_along(mb), mass_pp = mb, se_cluster_pp = sqrt(diag(V)) * 100, H_pp = 0.5 * sum(abs(mb))) }))
  list(res = res, strata = strata, masses = masses, Gs = Gs, keep = keep)
}
# --- gate: reproduce the stored reconstruction battery at K = 50 ---
rec50 <- moments(50L, dtm_ev, dtm_ev, 0L, "reconstruction")
stored <- as.data.table(x$tests_ho)[K == 50, .(test, stat, pval)]
cmp <- merge(rec50$res[, .(test, W_iid)], stored, by = "test")
cat("\nGATE reconstruction K=50 (recomputed vs stored):\n"); print(cmp)
stopifnot(all(abs(cmp$W_iid - cmp$stat) < 1e-6))
cat("GATE PASSED: fits, instruments and fold-in reproduce the stored statistics to 1e-6.\n")
# --- completion residuals: scoring half, mixture inferred from the fold-in half ---
spl <- split_tokens_binomial(dtm_ev, 0.5, seed = SEED_BASE)
cat("split: foldin/score objects:", paste(names(spl), collapse = ","), "; docs with positive scoring length:", sum(Matrix::rowSums(spl$score) > 0), "\n")
out <- list(rec50 = rec50)
for (k in K_test) out[[paste0("com", k)]] <- moments(k, spl$foldin, spl$score, 1L, "completion")
out$rec180 <- moments(180L, dtm_ev, dtm_ev, 0L, "reconstruction")
res <- rbindlist(lapply(out, `[[`, "res")); strata <- rbindlist(lapply(out, `[[`, "strata")); masses <- rbindlist(lapply(out, `[[`, "masses"))
cat("\n=== Moment tests: reconstruction vs completion residuals ===\n"); print(res[order(K, test, protocol)], digits = 4)
cat("\n=== Strata (Test 2 and Test 3 contrasts; Test 1) ===\n"); print(strata[order(K, test, protocol)], digits = 3)
cat("\n=== Frequency-group masses (Test 2 strata), pp of mean document mass, and H ===\n"); print(masses[test == "T2_freq_strata"][order(K, protocol, group)], digits = 3)
fwrite(res, file.path(SP, "stage3_mdna_completion_tests.csv")); fwrite(strata, file.path(SP, "stage3_mdna_completion_strata.csv")); fwrite(masses, file.path(SP, "stage3_mdna_completion_masses.csv"))
qs_save(lapply(out, function(o) list(Gs = o$Gs, keep = o$keep)), p_data("MDNA", "stage3_mdna_completion_G.qs2"))
cat(sprintf("done in %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
