# =============================================================================
# rescore_revision.R -- SHORT evaluation-only rescoring for the Sept-2026
# revision. Never estimates a model: options(optop.no_fit = TRUE) turns any
# fit-cache miss into an error (see prefit_pool() in R/source_all.R).
#
# Stages (run all by default, or pick with stages=tests,restarts,lemma,e1pilot):
#   tests    MD&A moment tests with K-SPECIFIC Test-3 strata at the reference
#            fit K_ref (+/- one grid step) and at the primary total-gain
#            selection K_sel (+/- one grid step). Moment matrices are saved, so
#            firm-cluster Wald tests come for free. Gates, in this order:
#              G1 forcing the pre-fix K = min(grid) strata must REPRODUCE the
#                 cached (buggy) Test-3 statistics  -> the recomputation path,
#                 the deterministic EM fold-in and the OpTop version in use
#                 all agree with the run that produced the paper's numbers;
#              G2 Tests 1-2 must reproduce the cached statistics (they were
#                 never affected by the defect);
#              G3 the mean residual vector at K_ref must reproduce the cached
#                 signed residual vocabulary.
#            Only then are the corrected Test-3 numbers written.
#   restarts the three production restarts at K_ref scored on ONE support
#            harmonised over the three fits (run_mdna_restarts.R scored each
#            on its own single-model support, which mixes optimisation
#            variability with a changing scoring construction).
#   lemma    Lemma S1 across code paths on the E6 seed-1 training fit, with the
#            count of zero-probability observed cells (unsmoothed WarpLDA phi).
#   wordfloor  how much the MD&A WORD-LEVEL indices depend on the expected-count
#            floor applied to zero-probability cells (document-level indices are
#            immune: such words always sit in the min-bin).
#   e1pilot  timing pilot for the deferred E1 rescoring (seed 1, Deviance,
#            reconstruction + completion), gated on reproducing doc_rep1.
#
# Each stage is timed; a stage is skipped when the elapsed budget is exhausted
# (pilot-first rule: nothing here may turn into a long run unnoticed).
#
# Usage:  Rscript Code/rescore_revision.R [stages=all] [K_ref=50] [K_sel=180]
#                                         [budget_min=15] [out_suffix=_rev1]
# Output: Data/MDNA/mdna_rescore<out_suffix>.qs2, Results/csv/*<out_suffix>.csv
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))
options(optop.no_fit = TRUE)
future::plan(future::sequential)

.args <- commandArgs(trailingOnly = TRUE)
P <- list(stages = "all", K_ref = 50L, K_sel = 180L, budget_min = 15,
          out_suffix = "_rev1")
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ", paste(names(P), collapse = ", "),
         call. = FALSE)
  P[[kv[1L]]] <- if (kv[1L] %in% c("stages", "out_suffix")) kv[2L]
                 else as.numeric(kv[2L])
}
STAGES <- if (P$stages == "all") c("tests", "restarts", "lemma", "wordfloor",
                                   "e1pilot") else
  strsplit(P$stages, ",", fixed = TRUE)[[1L]]
T0 <- proc.time()[["elapsed"]]
elapsed_min <- function() (proc.time()[["elapsed"]] - T0) / 60
in_budget <- function(what) {
  ok <- elapsed_min() < P$budget_min
  if (!ok) log_msg("BUDGET: %.1f min elapsed >= %.0f -- skipping %s",
                   elapsed_min(), P$budget_min, what)
  ok
}
gate <- function(ok, fmt, ...) {
  msg <- sprintf(fmt, ...)
  if (!isTRUE(ok)) stop("GATE FAILED: ", msg, call. = FALSE)
  log_msg("GATE ok: %s", msg)
}
rel_err <- function(a, b) max(abs(a - b) / pmax(abs(b), .Machine$double.xmin))
SFX <- P$out_suffix
SEED_BASE <- 1970L; C_PART <- 1
B_STRATA <- 5L; S_STRATA <- 5L; MIN_DOCFREQ <- 5L
out <- list(params = P)
f_out <- p_data("MDNA", sprintf("mdna_rescore%s.qs2", SFX))
# Saved after EVERY stage, merging with what earlier runs left on disk, so a
# failed gate in a later stage never discards a finished one.
save_out <- function() {
  if (file.exists(f_out)) {
    prev <- qs_read(f_out)
    for (nm in setdiff(names(prev), names(out))) out[[nm]] <<- prev[[nm]]
  }
  cache_put(out, f_out, list(experiment = "MDNA", profile = "revision"))
  log_msg("saved %s [%s]", basename(f_out), paste(setdiff(names(out), "params"),
                                                  collapse = ", "))
}

# ------------------------------ MD&A inputs -----------------------------------
need_mdna <- any(c("tests", "restarts", "wordfloor") %in% STAGES)
if (need_mdna) {
  prep <- qs_read(p_data("MDNA", "mdna_prep_2015_2016.qs2"))
  dtm_tr <- prep$dtm_train; dtm_ev <- prep$dtm_ev
  f_x0 <- p_data("MDNA", "mdna_results_MDNA_2015_2016.qs2")
  x0 <- qs_read(f_x0)                                                 # paper cache
  # this script quantifies the LEGACY word-level null against the paper's
  # definition, so its input must be the frozen pre-revision object (content hash)
  stopifnot(identical(word_null_convention_of(x0, f_x0),
                      WORD_NULL_CONVENTIONS[["legacy"]]))
  # identical corpus signature to run_mdna.R -> the production fit-cache entries
  sig <- list(corpus = "mdna_item7_pooled", y1 = 2015L, y2 = 2016L,
              dtm_hash = cfg_hash(list(dim(dtm_tr), Matrix::rowSums(dtm_tr)[1:20],
                                       colnames(dtm_tr)[1:50])))
  K_all <- x0$K_all
  cik <- data.table(doc_id = prep$dv_ev$doc_id,
                    cluster_id = as.character(prep$dv_ev$cik))
  stopifnot(identical(cik$doc_id, rownames(dtm_ev)))
}

# Cluster-robust Wald test of H0: E g = 0 from a J x q moment matrix: the
# variance of the mean is V = G/(G-1) * sum_g s_g s_g' / J^2 with s_g the
# within-cluster sum of centred moments (document-weighted mean preserved).
wald_cluster <- function(G, cl, label) {
  G <- as.matrix(G); J <- nrow(G); q <- ncol(G)
  gbar <- colMeans(G)
  S <- rowsum(sweep(G, 2L, gbar), cl)
  ng <- nrow(S)
  V <- ng / (ng - 1) * crossprod(S) / J^2
  W <- as.numeric(t(gbar) %*% solve(V) %*% gbar)
  data.table(test = label, stat_cluster = W, df = q,
             pval_cluster = pchisq(W, q, lower.tail = FALSE),
             n_clusters = ng, t_cluster_absmax = max(abs(gbar / sqrt(diag(V)))))
}

# Total mean residual mass by stratum (percentage points). Unlike the per-word
# contrasts, these are not rescaled by stratum size, so they measure how much
# probability is displaced BETWEEN groups (cancellation within a group and
# across documents is still possible). `strata` may contain NA = excluded words.
mass_by_stratum <- function(resid_mean, strata, label) {
  dt <- data.table(stratum = ifelse(is.na(strata), 0L, strata), r = resid_mean)
  o <- dt[, .(words = .N, mass_pp = 100 * sum(r)), by = stratum][order(stratum)]
  o[, `:=`(partition = label,
           stratum_lab = ifelse(stratum == 0L, "excluded", paste0("g", stratum)))]
  o[]
}

# ================================ stage: tests =================================
if ("tests" %in% STAGES && in_budget("tests")) {
  t_st <- proc.time()[["elapsed"]]
  step <- unique(diff(K_all))[1L]
  K_tests <- sort(unique(c(P$K_ref + c(-step, 0, step), P$K_sel + c(-step, 0, step))))
  K_tests <- as.integer(intersect(K_tests, K_all))
  K_bug <- min(K_all)                       # what the defect handed to every K
  log_msg("tests: K in {%s}; reference %d, primary selection %d",
          paste(K_tests, collapse = ", "), P$K_ref, P$K_sel)

  fits <- get_fits_cached(dtm_tr, K_all, "WarpLDA", 3L, SEED_BASE, sig)$models
  # training word scores exactly as production (all 20 fits in the partition)
  ins <- score_insample(fits, dtm_tr, C_PART, "dev",
                        word_at = sort(unique(c(K_bug, K_tests))))
  word_tr <- ins$word
  log_msg("tests: in-sample word scores done (%.1f min)", elapsed_min())
  Zs_by_K <- make_instruments_by_K(dtm_tr, word_tr, K_tests, B = B_STRATA,
                                   S = S_STRATA, min_docfreq = MIN_DOCFREQ)
  Z_bug <- make_instruments_by_K(dtm_tr, word_tr, K_bug, B = B_STRATA,
                                 S = S_STRATA,
                                 min_docfreq = MIN_DOCFREQ)[[1L]]$T3_fit_strata
  changed <- vapply(K_tests, function(k_v)
    sum(colSums(abs(Zs_by_K[[as.character(k_v)]]$T3_fit_strata - Z_bug)) > 0), 0)
  log_msg("tests: Test-3 instrument columns that differ from the K=%d strata: %s",
          K_bug, paste(sprintf("K=%d:%d", K_tests, as.integer(changed)),
                       collapse = "  "))

  blocks <- split(seq_len(nrow(dtm_ev)), ceiling(seq_len(nrow(dtm_ev)) / 1500L))
  tests <- list(); strata <- list(); moments <- list(); tests_cl <- list()
  resid_mean <- list(); timing <- list(); done_K <- integer(0)
  # essential fix first, then the new primary selection
  K_order <- c(intersect(c(P$K_ref, P$K_ref - step, P$K_ref + step), K_tests),
               intersect(c(P$K_sel, P$K_sel - step, P$K_sel + step), K_tests))
  for (k_v in as.integer(unique(K_order))) {
    if (!in_budget(sprintf("tests at K = %d", k_v))) break
    t_k <- proc.time()[["elapsed"]]
    k <- as.character(k_v)
    phi_k <- phi_from_fit(fits[[k]])
    Zs <- c(Zs_by_K[[k]], list(T3_bug_minK_strata = Z_bug))
    Gs <- lapply(Zs, function(Z) NULL); rs <- 0
    for (bi in seq_along(blocks)) {
      rows <- blocks[[bi]]
      th <- foldin_theta(fits[[k]], dtm_ev[rows, , drop = FALSE],
                         seed = SEED_BASE + bi)
      E <- resid_heldout(th, phi_k, dtm_ev[rows, , drop = FALSE])
      for (nm in names(Zs)) Gs[[nm]] <- rbind(Gs[[nm]], E %*% t(Zs[[nm]]))
      rs <- rs + colSums(E)
      rm(th, E); invisible(gc(FALSE))
    }
    mt <- lapply(names(Zs), function(nm) moment_test_from_G(Gs[[nm]], nm))
    tests[[k]] <- rbindlist(lapply(mt, `[[`, "result"), fill = TRUE)[, K := k_v]
    strata[[k]] <- rbindlist(lapply(mt, `[[`, "strata"))[, K := k_v]
    moments[[k]] <- lapply(Gs, function(G) { rownames(G) <- rownames(dtm_ev); G })
    tests_cl[[k]] <- rbindlist(lapply(names(Zs), function(nm)
      wald_cluster(Gs[[nm]], cik$cluster_id, nm)))[, K := k_v]
    resid_mean[[k]] <- rs / nrow(dtm_ev)
    done_K <- c(done_K, k_v)
    timing[[k]] <- data.table(K = k_v, sec = proc.time()[["elapsed"]] - t_k)
    log_msg("tests: K = %d done in %.0fs (%.1f min elapsed)", k_v,
            timing[[k]]$sec, elapsed_min())
  }
  tests <- rbindlist(tests, fill = TRUE); strata <- rbindlist(strata)
  tests_cl <- rbindlist(tests_cl); timing <- rbindlist(timing)

  # ---- gates against the paper cache (K_ref and its neighbours) ---------------
  old <- x0$tests_ho
  for (k_v in intersect(done_K, unique(old$K))) {
    o3 <- old[K == k_v & test == "T3_fit_strata"]
    n3 <- tests[K == k_v & test == "T3_bug_minK_strata"]
    gate(rel_err(n3$stat, o3$stat) < 1e-8 &&
           rel_err(n3$gbar_absmax, o3$gbar_absmax) < 1e-8,
         "G1 K=%d: forcing the K=%d strata reproduces the cached Test 3 (W = %.4f vs %.4f)",
         k_v, K_bug, n3$stat, o3$stat)
    for (tn in c("T1_freq_contrast", "T2_freq_strata")) {
      gate(rel_err(tests[K == k_v & test == tn, stat],
                   old[K == k_v & test == tn, stat]) < 1e-8,
           "G2 K=%d: %s reproduces the cached statistic (W = %.4f)", k_v, tn,
           old[K == k_v & test == tn, stat])
    }
  }
  if (P$K_ref %in% done_K) {
    rm_old <- x0$resid_words[match(colnames(dtm_ev), word), resid_mean]
    gate(max(abs(resid_mean[[as.character(P$K_ref)]] - rm_old)) < 1e-12,
         "G3 K=%d: mean residual vector reproduces the cached vocabulary (max |diff| %.2e)",
         P$K_ref, max(abs(resid_mean[[as.character(P$K_ref)]] - rm_old)))
    s1 <- strata[K == P$K_ref & test == "T3_fit_strata"][1L]
    log_msg("corrected Test 3 at K=%d: first contrast %.6e (cached, K=%d strata: %.6e)",
            P$K_ref, s1$gbar, K_bug,
            x0$strata_ho[K == P$K_ref & test == "T3_fit_strata"][1L]$gbar)
  }

  # ---- residual mass between vocabulary groups --------------------------------
  # Stratum membership is read back from the instrument matrices themselves
  # (row b is positive on stratum b, every row is negative on the reference
  # stratum, excluded words are zero), so it cannot drift from what was tested.
  strata_from_Z <- function(Z) {
    S_ <- nrow(Z) + 1L; s <- rep(NA_integer_, ncol(Z))
    for (b in seq_len(S_ - 1L)) s[Z[b, ] > 0] <- b
    s[Z[1L, ] < 0] <- S_
    s
  }
  stopifnot(identical(strata_from_Z(Zs_by_K[[1L]]$T2_freq_strata),
                      .freq_strata(dtm_tr, B_STRATA)))
  mass <- rbindlist(lapply(done_K, function(k_v) {
    k <- as.character(k_v)
    rbindlist(list(
      mass_by_stratum(resid_mean[[k]], strata_from_Z(Zs_by_K[[k]]$T2_freq_strata),
                      "training frequency (Test 2)"),
      mass_by_stratum(resid_mean[[k]], strata_from_Z(Zs_by_K[[k]]$T3_fit_strata),
                      "training word fit (Test 3)")
    ))[, K := k_v]
  }))
  mass[, half_abs_sum_pp := sum(abs(mass_pp)) / 2, by = .(K, partition)]

  resid_words <- rbindlist(lapply(done_K, function(k_v) data.table(
    K = k_v, word = colnames(dtm_ev), resid_mean = resid_mean[[as.character(k_v)]])))
  out$tests <- list(K_tests = done_K, K_bug = K_bug, tests = tests, strata = strata,
                    tests_cluster = tests_cl, moments = moments, mass = mass,
                    resid_words = resid_words, n_cols_changed = data.table(
                      K = K_tests, cols_changed = as.integer(changed)),
                    timing = timing)
  fwrite(tests, p_results("csv", sprintf("mdna_moment_tests%s.csv", SFX)))
  fwrite(strata, p_results("csv", sprintf("mdna_moment_strata%s.csv", SFX)))
  fwrite(tests_cl, p_results("csv", sprintf("mdna_moment_tests_cluster%s.csv", SFX)))
  fwrite(mass, p_results("csv", sprintf("mdna_residual_mass%s.csv", SFX)))
  log_msg("tests stage done in %.1f min", (proc.time()[["elapsed"]] - t_st) / 60)
  save_out()
  rm(fits, ins); invisible(gc())
}

# =============================== stage: restarts ===============================
if ("restarts" %in% STAGES && in_budget("restarts")) {
  t_st <- proc.time()[["elapsed"]]
  K_r <- as.integer(P$K_ref); N_STARTS_PROD <- 3L
  rf <- lapply(seq_len(N_STARTS_PROD), function(s)
    get_fits_cached(dtm_tr, K_r, "WarpLDA", 1L, 1969L + s, sig))
  prod_fit <- get_fits_cached(dtm_tr, K_r, "WarpLDA", 3L, SEED_BASE, sig)$models[[1L]]
  same <- vapply(rf, function(o) isTRUE(all.equal(
    as.matrix(o$models[[1L]]$tww), as.matrix(prod_fit$tww), tolerance = 1e-12)), NA)
  # SOFT check (WarpLDA need not be bit-reproducible across sessions): does one
  # single-start refit coincide with the production best-of-three fit?
  px_r <- vapply(rf, function(o) o$diagnostics$perplexity[1L], 0)
  log_msg(paste("restarts: single-start training perplexities %s; production",
                "best-of-3 %.4f; refits identical to the production fit: %s"),
          paste(sprintf("%.4f", px_r), collapse = " / "),
          x0$diagnostics[K == K_r, perplexity][1L],
          if (any(same)) paste(which(same), collapse = ",") else "NONE")
  if (sum(same) != 1L)
    log_msg("WARNING restarts: %d (not 1) refits equal the production fit",
            sum(same))
  # pseudo-K labels 1..3 so the three fits enter ONE harmonised partition
  fl <- lapply(rf, function(o) o$models[[1L]]); names(fl) <- as.character(1:3)
  pi_tr <- OpTop::optop_make_baseline(dtm_tr)$pi_glob
  rec <- score_heldout(fl, dtm_ev, dtm_ev, pi_tr, C_PART, "dev",
                       foldin_seed = SEED_BASE)
  spl <- split_tokens_binomial(dtm_ev, 0.5, seed = SEED_BASE)
  com <- score_heldout(fl, spl$foldin, spl$score, pi_tr, C_PART, "dev",
                       foldin_seed = SEED_BASE + 1L)
  disp <- rbindlist(list(
    data.table(protocol = "rec", rec$summary[, .(start = K, r2_macro, r2_micro, J_pos)]),
    data.table(protocol = "com", com$summary[, .(start = K, r2_macro, r2_micro, J_pos)])))
  disp[, `:=`(K = K_r, is_production = same[start],
              train_perplexity = vapply(rf, function(o) o$diagnostics$perplexity[1L],
                                        0)[start])]
  rng <- disp[, .(min = min(r2_macro), max = max(r2_macro),
                  range = max(r2_macro) - min(r2_macro)), by = protocol]
  # paired per-document comparison on the common support: is the dispersion
  # larger than evaluation noise?
  pr <- rbindlist(lapply(list(rec = rec, com = com), function(sc)
    paired_gains_all(sc$doc, family = "all_pairs")[, .(start_a = K, start_b = K_to,
                                                      delta_mean, se)]),
    idcol = "protocol")
  out$restarts <- list(dispersion = disp, range = rng, paired = pr,
                       support = "harmonised over the three restart fits only")
  fwrite(disp, p_results("csv", sprintf("mdna_restart_dispersion_common%s.csv", SFX)))
  print(disp); print(rng); print(pr)
  log_msg("restarts stage done in %.1f min", (proc.time()[["elapsed"]] - t_st) / 60)
  save_out()
  rm(rf, fl, rec, com, spl); invisible(gc())
}

# ================================ stage: lemma =================================
if ("lemma" %in% STAGES && in_budget("lemma")) {
  t_st <- proc.time()[["elapsed"]]
  f6 <- list.files(proj_path("Data", "E6"), "^e6_results_E6_full_.*warplda.*\\.qs2$",
                   full.names = TRUE)
  stopifnot(length(f6) >= 1L)
  e6 <- qs_read(f6[which.max(file.mtime(f6))])
  cfg6 <- e6$config
  rows <- list()
  for (sc_id in names(cfg6$scenarios)) {
    seeds <- make_seeds(cfg6$seed_base, 1L,
                        scenario_id = match(sc_id, names(cfg6$scenarios)))
    sc <- cfg6$scenarios[[sc_id]]
    set.seed(seeds$dgp_seed + 2L)
    Phi <- rdirichlet_mat(cfg6$K_true, cfg6$W, cfg6$beta_DGP)
    contam <- NULL
    if (isTRUE(sc$contam))
      contam <- list(pi_stop = make_stopword_dist(cfg6$W, cfg6$n_stop,
                                                  seed = seeds$dgp_seed + 71L),
                     w_mean = sc$w_mean, w_conc = 10, mode = "shared")
    sim_tr <- sim_lda_corpus(cfg6$J_train, cfg6$W, cfg6$K_true, cfg6$alpha_DGP,
                             cfg6$beta_DGP, cfg6$length_spec, seed = seeds$dgp_seed,
                             Phi = Phi, contamination = contam, doc_prefix = "tr")
    fit <- get_fits_cached(sim_tr$dtm, cfg6$K_true, cfg6$fit_method, cfg6$n_starts,
                           seeds$fit_seed_base,
                           dgp_signature(cfg6, seeds, scenario = sc_id))$models
    w <- score_insample(fit, sim_tr$dtm, cfg6$c, "dev", word_at = cfg6$K_true)$word
    w <- w[match(colnames(sim_tr$dtm), word_id)]
    l2 <- lemma2_residual(theta_from_fit(fit[[1L]]), phi_from_fit(fit[[1L]]),
                          sim_tr$dtm, word_d_model = w$d_model)
    rows[[sc_id]] <- data.table(scenario = sc_id, dev_doc = l2$dev_doc,
                                dev_word_pkg = l2$dev_word_pkg,
                                rel_resid_pkg = l2$rel_resid_pkg,
                                max_word_diff = l2$max_word_diff,
                                reorder_resid = l2$resid,
                                phi_zero_share = mean(phi_from_fit(fit[[1L]]) == 0),
                                n_cells = l2$n_cells, n_floored = l2$n_floored,
                                token_share_floored = l2$token_share_floored,
                                dev_share_floored = l2$dev_share_floored,
                                n_words_floored = l2$n_words_floored)
  }
  lem <- rbindlist(rows)
  print(lem)
  gate(all(abs(lem$rel_resid_pkg) < 1e-8),
       "lemma: doc-wise total == package word-wise total (max rel %.2e)",
       max(abs(lem$rel_resid_pkg)))
  out$lemma <- lem
  fwrite(lem, p_results("csv", sprintf("e6_lemma_crosspath%s.csv", SFX)))
  print(lem)
  log_msg("lemma stage done in %.1f min", (proc.time()[["elapsed"]] - t_st) / 60)
  save_out()
}

# ============================== stage: wordfloor ===============================
if ("wordfloor" %in% STAGES && in_budget("wordfloor")) {
  t_st <- proc.time()[["elapsed"]]
  K_w <- as.integer(P$K_ref)
  fitw <- get_fits_cached(dtm_tr, K_w, "WarpLDA", 3L, SEED_BASE, sig)$models[[1L]]
  phw <- phi_from_fit(fitw); tphw <- t(phw)
  pi_tr <- as.numeric(OpTop::optop_make_baseline(dtm_tr)$pi_glob)
  Lw <- as.numeric(Matrix::rowSums(dtm_ev))
  blocks <- split(seq_len(nrow(dtm_ev)), ceiling(seq_len(nrow(dtm_ev)) / 1500L))
  cells <- list(); Ltheta <- numeric(K_w)
  for (bi in seq_along(blocks)) {
    rows <- blocks[[bi]]
    th <- foldin_theta(fitw, dtm_ev[rows, , drop = FALSE], seed = SEED_BASE)
    Ltheta <- Ltheta + as.numeric(crossprod(th, Lw[rows]))
    Tb <- as(as(dtm_ev[rows, , drop = FALSE], "generalMatrix"), "TsparseMatrix")
    ib <- Tb@i + 1L; jb <- Tb@j + 1L
    pz <- numeric(length(ib))
    for (ch in split(seq_along(ib), ceiling(seq_along(ib) / 250000L)))
      pz[ch] <- rowSums(th[ib[ch], , drop = FALSE] * tphw[jb[ch], , drop = FALSE])
    cells[[bi]] <- data.table(j = jb, N = Tb@x, E = pz * Lw[rows][ib],
                              B = pi_tr[jb] * Lw[rows][ib])
    rm(th, Tb); invisible(gc(FALSE))
  }
  cells <- rbindlist(cells)
  E_w <- as.numeric(Ltheta %*% phw); N_w <- as.numeric(Matrix::colSums(dtm_ev))
  B_w <- pi_tr * sum(Lw)
  word_dev <- function(e_floor, drop = FALSE) {
    cc <- if (drop) cells[E >= 1e-12] else cells
    a <- cc[, .(s = sum(2 * N * log(N / pmax(E, e_floor)))), by = j]
    d <- numeric(ncol(dtm_ev)); d[a$j] <- a$s
    d - 2 * (N_w - E_w)
  }
  # --- null deviance: the package convention vs the paper's definition ---------
  # OpTop (<= 0.20.1, src/index_core.cpp; 0.20.1 adds a validity warning, the
  # kernel is unchanged) returns the word-level FITTED deviance
  # in Poisson form, 2 sum_j [N log(N/E) - (N - E)], but the word-level NULL
  # deviance WITHOUT its linear term, 2 sum_j N log(N/B). In-sample the omitted
  # term is identically zero (B_w = N_w because the baseline is the corpus
  # marginal), so nothing in-sample is affected -- including the training word
  # scores behind Test 3. HELD-OUT, with a TRAINING baseline and EVALUATION
  # counts, N_w - B_w != 0 and the two conventions differ. Section 2.4 of the
  # paper defines the null deviance in Poisson form, so both are reported.
  a0 <- cells[, .(s = sum(2 * N * log(N / B))), by = j]
  dn_pkg <- numeric(ncol(dtm_ev)); dn_pkg[a0$j] <- a0$s
  dn_paper <- dn_pkg - 2 * (N_w - B_w)
  ws <- x0$wstar[match(colnames(dtm_ev), word_id)]
  stopifnot(all(ws$K == K_w))
  d12 <- word_dev(1e-12)
  gate(rel_err(d12[ws$keep], ws$d_model[ws$keep]) < 1e-8,
       "wordfloor A: fitted word deviances (E floored at 1e-12, Poisson form) reproduce the cached package values (max rel %.1e)",
       rel_err(d12[ws$keep], ws$d_model[ws$keep]))
  gate(rel_err(dn_pkg[ws$keep], ws$d_null[ws$keep]) < 1e-8,
       "wordfloor B: cached word-level NULL deviances equal 2 sum N log(N/B) WITHOUT the Poisson linear term (max rel %.1e)",
       rel_err(dn_pkg[ws$keep], ws$d_null[ws$keep]))
  zc <- cells[E < 1e-12, .(n_zero_cells = .N, tokens = sum(N)), by = j]
  has_zero <- seq_len(ncol(dtm_ev)) %in% zc$j
  elig <- ws$doc_freq >= MIN_DOCFREQ & ws$B_w >= 5     # the Section-2.4 filter
  summ <- function(dm, dn, label_floor, label_null) {
    k <- elig & dn > 0
    r2 <- 1 - dm / dn
    data.table(null_deviance = label_null, fitted_floor = label_floor, n_words = sum(k),
               w_micro = 1 - sum(dm[k]) / sum(dn[k]), w_macro = mean(r2[k]),
               gap = (1 - sum(dm[k]) / sum(dn[k])) - mean(r2[k]),
               w_macro_excl_zero_cell_words = mean(r2[k & !has_zero]),
               min_r2 = min(r2[k]), n_r2_below_m1 = sum(r2[k] < -1))
  }
  floors <- list("E floor 1e-8" = word_dev(1e-8), "E floor 1e-12 (package)" = d12,
                 "E floor 1e-16" = word_dev(1e-16),
                 "zero-probability cells dropped" = word_dev(1e-12, drop = TRUE))
  fl <- rbindlist(c(
    lapply(names(floors), function(nm) summ(floors[[nm]], dn_pkg, nm,
                                            "package: no linear term")),
    lapply(names(floors), function(nm) summ(floors[[nm]], dn_paper, nm,
                                            "paper: Poisson form"))))
  cw <- x0$wcurve[K == K_w]
  base_row <- fl[null_deviance == "package: no linear term" &
                   fitted_floor == "E floor 1e-12 (package)"]
  gate(abs(base_row$w_micro - cw$r2_micro_word) < 1e-8 &&
         abs(base_row$w_macro - cw$r2_macro_word) < 1e-8 && base_row$n_words == cw$n_words,
       "wordfloor C: package conventions reproduce the cached word-level Micro %.4f / Macro %.4f (%d words) at K=%d",
       cw$r2_micro_word, cw$r2_macro_word, cw$n_words, K_w)
  # Micro over the WHOLE grid under the paper's null: sum_w D_w(K) is recovered
  # from the cached Micro and the package null total on the same word set.
  k0 <- elig & dn_pkg > 0
  micro_all <- x0$wcurve[, .(K, w_micro_package = r2_micro_word,
    w_micro_paper_null = 1 - (1 - r2_micro_word) * sum(dn_pkg[k0]) / sum(dn_paper[k0]),
    w_macro_package = r2_macro_word)]
  worst <- data.table(word = colnames(dtm_ev), doc_freq = ws$doc_freq,
                      r2_package = 1 - d12 / dn_pkg, r2_paper_null = 1 - d12 / dn_paper,
                      r2_paper_null_dropzero = 1 - word_dev(1e-12, drop = TRUE) / dn_paper,
                      zero_cells = 0L, elig = elig)
  worst[zc$j, zero_cells := zc$n_zero_cells]
  worst_pkg <- worst[elig & dn_pkg > 0][order(r2_package)][1:40]
  worst_pap <- worst[elig & dn_paper > 0][order(r2_paper_null)][1:40]
  diag <- data.table(K = K_w, n_cells = nrow(cells), n_zero_cells = cells[, sum(E < 1e-12)],
                     token_share_zero = cells[E < 1e-12, sum(N)] / sum(cells$N),
                     n_words_with_zero = sum(has_zero),
                     n_elig_words_with_zero = sum(has_zero & elig),
                     phi_zero_share = mean(phw == 0),
                     worst20_pkg_with_zero = worst_pkg[1:20, sum(zero_cells > 0)],
                     worst20_paper_with_zero = worst_pap[1:20, sum(zero_cells > 0)],
                     worst20_overlap = length(intersect(worst_pkg$word[1:20], worst_pap$word[1:20])),
                     rel_change_null_median = median(abs(dn_paper[k0] / dn_pkg[k0] - 1)),
                     n_elig_pkg = sum(k0), n_elig_paper = sum(elig & dn_paper > 0))
  out$wordfloor <- list(diag = diag, conventions = fl, micro_over_K = micro_all,
                        worst_package = worst_pkg, worst_paper = worst_pap)
  fwrite(fl, p_results("csv", sprintf("mdna_wordlevel_conventions%s.csv", SFX)))
  fwrite(micro_all, p_results("csv", sprintf("mdna_wordlevel_micro_over_K%s.csv", SFX)))
  fwrite(worst_pap, p_results("csv", sprintf("mdna_wordlevel_worst_fit_paper_null%s.csv", SFX)))
  fwrite(worst_pkg, p_results("csv", sprintf("mdna_wordlevel_worst_fit_package%s.csv", SFX)))
  print(diag); print(fl, digits = 4); print(micro_all[K %in% c(10, 50, 100, 150, 200)], digits = 4)
  cat("worst fit, paper null:\n"); print(worst_pap[1:10, .(word, doc_freq, r2_paper_null, r2_package, zero_cells)], digits = 4)
  cat("worst fit, package:\n"); print(worst_pkg[1:10, .(word, doc_freq, r2_package, r2_paper_null, zero_cells)], digits = 4)
  log_msg("wordfloor stage done in %.1f min", (proc.time()[["elapsed"]] - t_st) / 60)
  save_out()
  rm(cells); invisible(gc())
}

# =============================== stage: e1pilot ================================
if ("e1pilot" %in% STAGES && in_budget("e1pilot")) {
  t_st <- proc.time()[["elapsed"]]
  f1 <- p_data("E1", "e1_results_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2")
  e1 <- qs_read(f1); cfg1 <- e1$config
  seeds <- make_seeds(cfg1$seed_base, 1L)
  sim_tr <- sim_lda_corpus(cfg1$J_train, cfg1$W, cfg1$K_true, cfg1$alpha_DGP,
                           cfg1$beta_DGP, cfg1$length_spec, seed = seeds$dgp_seed,
                           doc_prefix = "tr")
  sim_ev <- sim_lda_corpus(cfg1$J_eval, cfg1$W, cfg1$K_true, cfg1$alpha_DGP,
                           cfg1$beta_DGP, cfg1$length_spec,
                           seed = seeds$dgp_seed + 7L, Phi = sim_tr$Phi,
                           doc_prefix = "ev")
  fits1 <- get_fits_cached(sim_tr$dtm, cfg1$K_grid, cfg1$fit_method, cfg1$n_starts,
                           seeds$fit_seed_base, dgp_signature(cfg1, seeds))$models
  pi1 <- OpTop::optop_make_baseline(sim_tr$dtm)$pi_glob
  t_r <- proc.time()[["elapsed"]]
  rec1 <- score_heldout(fits1, sim_ev$dtm, sim_ev$dtm, pi1, cfg1$c, "dev",
                        foldin_seed = seeds$split_seed)
  sec_rec <- proc.time()[["elapsed"]] - t_r; t_c <- proc.time()[["elapsed"]]
  spl1 <- split_tokens_binomial(sim_ev$dtm, cfg1$completion_prop,
                                seed = seeds$split_seed)
  com1 <- score_heldout(fits1, spl1$foldin, spl1$score, pi1, cfg1$c, "dev",
                        foldin_seed = seeds$split_seed + 1L)
  sec_com <- proc.time()[["elapsed"]] - t_c
  chk <- function(new, mode) {
    o <- e1$doc_rep1[eval == mode & metric == "dev", .(K, doc_id, r2_old = r2_doc)]
    m <- merge(new$doc[metric == "dev", .(K, doc_id, r2_doc)], o, by = c("K", "doc_id"))
    c(n = nrow(m), na_mismatch = sum(is.na(m$r2_doc) != is.na(m$r2_old)),
      max_abs = max(abs(m$r2_doc - m$r2_old), na.rm = TRUE))
  }
  c_rec <- chk(rec1, "ho_reconstruction"); c_com <- chk(com1, "ho_completion")
  gate(c_rec[["na_mismatch"]] == 0 && c_com[["na_mismatch"]] == 0 &&
         c_rec[["max_abs"]] < 1e-8 && c_com[["max_abs"]] < 1e-8,
       "e1pilot: seed-1 rescoring reproduces the cached document scores (max |diff| rec %.1e, com %.1e)",
       c_rec[["max_abs"]], c_com[["max_abs"]])
  sel1 <- rbindlist(list(
    select_k_all_rules(rec1$doc, cfg1$eps_grid, cfg1$sel_alpha)[, protocol := "rec"],
    select_k_all_rules(com1$doc, cfg1$eps_grid, cfg1$sel_alpha)[, protocol := "com"]))
  out$e1pilot <- list(sec_rec = sec_rec, sec_com = sec_com, selections = sel1,
                      check = rbind(rec = c_rec, com = c_com))
  log_msg(paste("e1pilot: ONE seed, Deviance only -- reconstruction %.0fs + completion",
                "%.0fs = %.1f min on one core; 10 seeds ~ %.0f core-min,",
                "100 selection replicates ~ %.0f core-min"),
          sec_rec, sec_com, (sec_rec + sec_com) / 60, 10 * (sec_rec + sec_com) / 60,
          100 * (sec_rec + sec_com) / 60)
  print(sel1[, .(protocol, rule, eps, K_hat, certified, n_comparisons, max_ub)])
  log_msg("e1pilot stage done in %.1f min", (proc.time()[["elapsed"]] - t_st) / 60)
  save_out()
}

log_msg("=== rescore_revision complete | total %.1f min ===", elapsed_min())
