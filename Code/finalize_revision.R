# =============================================================================
# finalize_revision.R -- acceptance checks on the outputs of the revision batch.
#
# Reads the EXPLICIT output paths the batch runner recorded
# (Results/revision<suffix>/inputs.json); nothing is located by pattern or by
# modification time. Every check is recorded, none stops the script: the verdict
# is the whole table (Results/revision<suffix>/acceptance.json and .csv) and the
# exit status is non-zero unless every check passes.
#
# What is accepted:
#   A  every production output exists, is stamped with the scoring convention
#      (Poisson-form held-out word null on the scored tokens) and the package;
#   B  its configuration is the source configuration, up to declared overrides;
#   C  the held-out baseline deviance of the word "condition" on the MD&A
#      evaluation corpus is 9562.168059 (legacy log-only form: 8771.558890),
#      recomputed here from the definition, independently of the pipeline;
#   D  independent fitted / null discrepancy checks (E6 Lemma S1 across code
#      paths; MD&A word summaries at K = 50 against the verified rev1 values);
#   E  scores the revision does NOT touch reproduce their historical values
#      (document-level indices, Tests 1-2, Test 3 at the smallest tested K,
#      pointwise selections), and the corrected MD&A Test 3 equals the verified
#      values at K in {40, 50, 60, 170, 180, 190};
#   F  the fit cache is byte-identical to the pre-revision manifest.
#
# Usage: Rscript Code/finalize_revision.R suffix=_rev2 [skip_fits=0]
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))

.args <- commandArgs(trailingOnly = TRUE)
arg <- function(key, default) {
  v <- sub(paste0("^", key, "="), "", grep(paste0("^", key, "="), .args, value = TRUE))
  if (length(v)) v[[1L]] else default
}
SFX <- arg("suffix", "_rev2"); SKIP_FITS <- arg("skip_fits", "0") == "1"
dest <- proj_path("Results", paste0("revision", SFX))
man_f <- file.path(dest, "inputs.json")
if (!file.exists(man_f)) stop("no ", man_f, ": run the batch first", call. = FALSE)
man <- jsonlite::fromJSON(man_f, simplifyVector = FALSE)
out_of <- function(key) {
  p <- unlist(man$outputs[[key]])
  if (length(p) != 1L) stop("stage '", key, "' declares no single output", call. = FALSE)
  proj_path(p)
}
TOL <- 1e-8
checks <- list()
record <- function(id, what, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(id = id, check = what,
                                               pass = isTRUE(ok), detail = detail)
  log_msg("%s  %-4s %s%s", if (isTRUE(ok)) "ok  " else "FAIL", id, what,
          if (nzchar(detail)) paste0("  [", detail, "]") else "")
}
try_check <- function(id, what, expr) {
  r <- tryCatch(expr, error = function(e) list(ok = FALSE, detail = conditionMessage(e)))
  record(id, what, r$ok, r$detail %||% "")
}
maxdiff <- function(a, b) if (length(a) == length(b) && length(a)) max(abs(a - b), na.rm = TRUE) else Inf
src_of <- function(path) sub(paste0(SFX, "\\.qs2$"), ".qs2", path)

sim_keys <- c("e1_base", "e1_dgp2", "e2_base", "e2_dgp2", "e6", "e4_base", "e4_prior")
mdna_keys <- c("mdna", "mdna_c05", "mdna_c2", "mdna_g100")
all_keys <- c(mdna_keys, sim_keys, "selreps", "e5_resolution", "unbinned", "restarts")

# ---- A. existence + stamps -------------------------------------------------------
objs <- list()
for (k in all_keys) {
  p <- out_of(k)
  record("A1", sprintf("output of '%s' exists", k), file.exists(p), basename(p))
  if (file.exists(p)) objs[[k]] <- qs_read(p)
}
for (k in names(objs)) {
  sc <- attr(objs[[k]], "run_meta")$scoring
  record("A2", sprintf("'%s' stamped: Poisson word null on scored tokens, OpTop >= 0.20.1", k),
         identical(sc$word_null_convention, "poisson_scored_tokens") &&
           utils::compareVersion(sc$optop_version %||% "0", "0.20.1") >= 0,
         sprintf("%s / OpTop %s", sc$word_null_convention %||% "unstamped",
                 sc$optop_version %||% "?"))
}

# ---- B. configuration ---------------------------------------------------------------
for (k in intersect(sim_keys, names(objs))) try_check(
  "B1", sprintf("'%s' configuration equals its source (declared overrides aside)", k), {
    c0 <- qs_read(src_of(out_of(k)))$config; c1 <- objs[[k]]$config
    drop <- c("J_truth", "K_test")      # J_truth: declared override; K_test: normalised by run_E4.R
    d <- all.equal(c0[setdiff(names(c0), drop)], c1[setdiff(names(c1), drop)])
    list(ok = isTRUE(d), detail = if (isTRUE(d)) "" else paste(head(d, 2L), collapse = "; "))
  })
for (k in intersect(c("e2_base", "e2_dgp2"), names(objs)))
  record("B2", sprintf("'%s' reference sample has 20,000 documents, with stored SEs", k),
         identical(as.numeric(objs[[k]]$config$J_truth), 20000) &&
           all(c("mu_se", "delta_true_se", "gap_true_se") %in% names(objs[[k]]$truth)),
         sprintf("J_truth = %s", objs[[k]]$config$J_truth))
if (!is.null(objs$mdna)) {
  ic <- objs$mdna$input_check
  for (k in intersect(mdna_keys, names(objs)))
    record("B3", sprintf("'%s': delta = 1, unrefined grid, same evaluation data", k),
           identical(as.numeric(objs[[k]]$input_check$min_null), 1) &&
             identical(objs[[k]]$input_check$dtm_ev_hash, ic$dtm_ev_hash) &&
             identical(as.integer(objs[[k]]$params$refine_span), 0L),
           sprintf("c = %s, |K| = %d", objs[[k]]$input_check$c_part, length(objs[[k]]$K_all)))
  record("B4", "MD&A diagnostics kept at K in {40,50,60,170,180,190}",
         all(c(40L, 50L, 60L, 170L, 180L, 190L) %in% objs$mdna$K_test),
         paste(objs$mdna$K_test, collapse = ","))
}

# ---- C. the "condition" acceptance value ---------------------------------------------
if (!is.null(objs$mdna)) try_check(
  "C1", "D_null('condition') = 9562.168059 from the definition (legacy 8771.558890)", {
    prep <- qs_read(p_data("MDNA", "mdna_prep_2015_2016.qs2"))
    n <- as.numeric(prep$dtm_ev[, "condition"]); L <- as.numeric(Matrix::rowSums(prep$dtm_ev))
    pi_c <- as.numeric(OpTop::optop_make_baseline(prep$dtm_train)$pi_glob[
      match("condition", colnames(prep$dtm_train))])
    b <- L * pi_c
    log_only <- 2 * sum(ifelse(n > 0, n * log(n / b), 0))
    poisson <- log_only - 2 * (sum(n) - sum(b))
    got <- objs$mdna$word_rec[word_id == "condition", d_null]
    list(ok = abs(poisson - 9562.168059) < 1e-5 && abs(log_only - 8771.558890) < 1e-5 &&
           length(got) == length(objs$mdna$K_all) && max(abs(got - poisson)) < 1e-6,
         detail = sprintf("definition %.6f | pipeline %.6f (all %d K) | log-only %.6f",
                          poisson, got[1L], length(got), log_only))
  })

# ---- D. independent discrepancy checks --------------------------------------------------
if (!is.null(objs$e6)) {
  l2 <- objs$e6$lemma2
  record("D1", "E6 Lemma S1: document total == package word totals, every unit",
         all(abs(l2$rel_resid_pkg) < 1e-10) && nrow(l2) == 2L * objs$e6$config$S,
         sprintf("max |rel| %.1e over %d units", max(abs(l2$rel_resid_pkg)), nrow(l2)))
  raw <- objs$e6$word_raw
  record("D2", "E6 raw word discrepancies kept for every K, both scenarios, both samples",
         !is.null(raw) && setequal(unique(raw$K), objs$e6$config$K_grid) &&
           uniqueN(raw$scenario) == length(objs$e6$config$scenarios) &&
           uniqueN(raw[["eval"]]) == 2L,
         sprintf("%s rows", format(nrow(raw), big.mark = ",")))
}
if (!is.null(objs$mdna)) try_check(
  "D3", "MD&A word-level summaries at K = 50 equal the verified rev1 values", {
    r1 <- qs_read(p_data("MDNA", "mdna_rescore_rev1.qs2"))$wordfloor$conventions
    ref <- r1[null_deviance != "package: no linear term"][1L]
    now <- objs$mdna$wcurve[K == 50L & metric == "dev"]
    list(ok = abs(now$r2_micro_word - ref$w_micro) < 1e-6 &&
           abs(now$r2_macro_word - ref$w_macro) < 1e-6 && now$n_words == ref$n_words,
         detail = sprintf("Micro %.4f / Macro %.4f on %d words", now$r2_micro_word,
                          now$r2_macro_word, now$n_words))
  })
if (!is.null(objs$mdna))
  record("D4", "MD&A raw word discrepancies kept at every K (held-out and training)",
         setequal(unique(objs$mdna$word_rec$K), objs$mdna$K_all) &&
           setequal(unique(objs$mdna$word_train$K), objs$mdna$K_all),
         sprintf("%s held-out rows", format(nrow(objs$mdna$word_rec), big.mark = ",")))

# ---- E. untouched scores reproduce; corrected Test 3 equals the verified values ---------
if (!is.null(objs$mdna)) {
  x0 <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016.qs2"))
  try_check("E1", "MD&A document-level indices reproduce the published cache", {
    m <- merge(objs$mdna$summary, x0$summary, by = c("K", "metric", "eval"),
               suffixes = c("", ".old"))
    list(ok = nrow(m) == nrow(x0$summary) &&
           max(abs(m$r2_macro - m$r2_macro.old), abs(m$r2_micro - m$r2_micro.old)) < TOL,
         detail = sprintf("max |diff| %.1e over %d rows",
                          max(abs(m$r2_macro - m$r2_macro.old),
                              abs(m$r2_micro - m$r2_micro.old)), nrow(m)))
  })
  try_check("E2", "MD&A Tests 1-2 reproduce the published cache (K = 40, 50, 60)", {
    m <- merge(objs$mdna$tests_ho[test != "T3_fit_strata"],
               x0$tests_ho[test != "T3_fit_strata"], by = c("K", "test"),
               suffixes = c("", ".old"))
    list(ok = nrow(m) == 6L && max(abs(m$stat / m$stat.old - 1)) < TOL,
         detail = sprintf("max rel diff %.1e", max(abs(m$stat / m$stat.old - 1))))
  })
  try_check("E3", "MD&A corrected Test 3 equals the verified values at the six K", {
    r1 <- qs_read(p_data("MDNA", "mdna_rescore_rev1.qs2"))$tests$tests
    m <- merge(objs$mdna$tests_ho[test == "T3_fit_strata"],
               r1[test == "T3_fit_strata"], by = c("K", "test"), suffixes = c("", ".v"))
    list(ok = nrow(m) == 6L && max(abs(m$stat / m$stat.v - 1)) < 1e-6,
         detail = paste(sprintf("K=%d: %.2f", m$K, m$stat), collapse = ", "))
  })
  try_check("E4", "MD&A selections under the three rules equal the verified values", {
    s <- objs$mdna$sel_all[metric == "dev" & eps == 0.01]
    g <- function(ev, rl) s[s[["eval"]] == ev & rule == rl, K_hat]
    got <- c(g("com", "adjacent_pointwise"), g("com", "adjacent_simultaneous"),
             g("com", "total_gain"), g("rec", "adjacent_pointwise"),
             g("rec", "adjacent_simultaneous"), g("rec", "total_gain"))
    list(ok = identical(as.integer(got), c(50L, 80L, 180L, 80L, 80L, 180L)),
         detail = paste(got, collapse = "/"))
  })
}
for (k in intersect(c("e1_base", "e1_dgp2"), names(objs))) try_check(
  "E5", sprintf("'%s' indices and pointwise selections reproduce the published cache", k), {
    o <- qs_read(src_of(out_of(k)))
    m <- merge(objs[[k]]$summary, o$summary, by = c("K", "metric", "eval", "replicate"),
               suffixes = c("", ".old"))
    kh <- merge(objs[[k]]$khat, o$khat, by = c("metric", "eval", "replicate", "rule", "eps"),
                suffixes = c("", ".old"))
    d <- max(abs(m$r2_macro - m$r2_macro.old), abs(m$r2_micro - m$r2_micro.old), na.rm = TRUE)
    list(ok = nrow(m) == nrow(o$summary) && d < TOL &&
           identical(kh$K_hat, kh$K_hat.old) &&
           uniqueN(objs[[k]]$sel_all$replicate) == o$config$S,
         detail = sprintf("max |diff| %.1e; %d replicates under all rules", d,
                          uniqueN(objs[[k]]$sel_all$replicate)))
  })
for (k in intersect(c("e2_base", "e2_dgp2"), names(objs))) try_check(
  "E6", sprintf("'%s' evaluation-replicate estimates reproduce the published cache", k), {
    o <- qs_read(src_of(out_of(k)))
    by_ <- c("K", "J_ev", "train_seed", "rep")
    m <- merge(objs[[k]]$cover, o$cover, by = by_, suffixes = c("", ".old"))
    d <- max(abs(m$r2_macro - m$r2_macro.old), abs(m$se - m$se.old))
    list(ok = nrow(m) == nrow(o$cover) && d < TOL,
         detail = sprintf("max |diff| %.1e over %d replicate rows (targets differ by design)",
                          d, nrow(m)))
  })
for (k in intersect(c("e4_base", "e4_prior"), names(objs))) try_check(
  "E7", sprintf("'%s': Tests 1-2 and Test 3 at the smallest K reproduce; Test 3 at K* changes", k), {
    o <- qs_read(src_of(out_of(k)))
    by_ <- c("test", "K", "rep", "train_seed")
    m <- merge(objs[[k]]$size, o$size, by = by_, suffixes = c("", ".old"))
    kmin <- min(m$K)
    same <- m[test != "T3_fit_strata" | K == kmin]
    chg <- m[test == "T3_fit_strata" & K != kmin]
    d <- max(abs(same$stat / same$stat.old - 1))
    list(ok = nrow(m) == nrow(o$size) && d < 1e-6 &&
           mean(abs(chg$stat / chg$stat.old - 1) > 1e-6) > 0.9 &&
           length(objs[[k]]$centers) == o$config$S_train && length(objs[[k]]$moments) > 0L,
         detail = sprintf("unchanged rows: max rel diff %.1e; Test 3 at K > %d changed in %.0f%% of %d rows",
                          d, kmin, 100 * mean(abs(chg$stat / chg$stat.old - 1) > 1e-6), nrow(chg)))
  })
if (!is.null(objs$e6)) try_check(
  "E8", "E6 in-sample word curves reproduce the published cache", {
    o <- qs_read(src_of(out_of("e6")))
    by_ <- c("K", "metric", "eval", "scenario", "replicate")
    a <- objs$e6$curve; a <- a[a[["eval"]] == "insample"]
    b <- o$curve; b <- b[b[["eval"]] == "insample"]
    m <- merge(a, b, by = by_, suffixes = c("", ".old"))
    d <- max(abs(m$r2_micro_word - m$r2_micro_word.old), abs(m$r2_macro_word - m$r2_macro_word.old))
    list(ok = nrow(m) == nrow(b) && d < TOL, detail = sprintf("max |diff| %.1e", d))
  })
if (!is.null(objs$selreps)) try_check(
  "E9", "selreps: pointwise selections reproduce; all three rules for every replication", {
    o <- qs_read(src_of(out_of("selreps")))
    n_rep <- objs$selreps$config$S * objs$selreps$R_sel
    sa <- objs$selreps$sel_all
    list(ok = isTRUE(all.equal(objs$selreps$khat, o$khat, check.attributes = FALSE)) &&
           all(sa[, .N, by = c("metric", "eval", "rule", "eps")]$N == n_rep),
         detail = sprintf("%d replications x %d rules", n_rep, uniqueN(sa$rule)))
  })
if (!is.null(objs$restarts))
  record("E10", "restarts scored on the support of the production grid and the restarts",
         identical(as.integer(objs$restarts$params$common), 2L) &&
           nrow(objs$restarts$dispersion) == 6L,
         objs$restarts$support %||% "")
if (!is.null(objs$unbinned))
  record("E11", "unbinned comparator: three floors, both protocols, every candidate",
         setequal(objs$unbinned$floors, c(1e-14, 1e-12, 1e-10)) &&
           all(objs$unbinned$summary[, .N, by = c("source", "K")]$N == 6L) &&
           all(is.finite(objs$unbinned$summary$log_score_per_token)),
         sprintf("%d sources", uniqueN(objs$unbinned$summary$source)))
for (k in intersect(c("mdna", "e1_base", "e5_resolution"), names(objs)))
  record("E12", sprintf("'%s' carries support-resolution rows", k),
         !is.null(objs[[k]]$resolution) && nrow(objs[[k]]$resolution) > 0L &&
           all(c("retained_cells", "observed_pooled_count", "predicted_pooled_share") %in%
                 names(objs[[k]]$resolution)),
         sprintf("%s rows", format(nrow(objs[[k]]$resolution), big.mark = ",")))

# ---- F. fit cache ---------------------------------------------------------------------------
if (!SKIP_FITS) try_check("F1", "every cached fit is byte-identical to the pre-revision manifest", {
  mf <- fread(p_results("csv", "fit_manifest_rev1.csv"))
  now <- vapply(proj_path(mf$path), function(p)
    if (file.exists(p)) digest::digest(p, algo = "sha256", file = TRUE) else NA_character_, "")
  list(ok = !anyNA(now) && all(now == mf$sha256),
       detail = sprintf("%d files; %d missing, %d changed", nrow(mf), sum(is.na(now)),
                        sum(now != mf$sha256, na.rm = TRUE)))
})

# ---- verdict -----------------------------------------------------------------------------------
tab <- rbindlist(checks)
fwrite(tab, file.path(dest, "acceptance.csv"))
jsonlite::write_json(list(suffix = SFX, at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                          n_checks = nrow(tab), n_failed = sum(!tab$pass),
                          accepted = all(tab$pass), checks = tab),
                     file.path(dest, "acceptance.json"), pretty = TRUE, auto_unbox = TRUE)
log_msg("=== acceptance: %d checks, %d failed -> %s ===", nrow(tab), sum(!tab$pass),
        if (all(tab$pass)) "ACCEPTED" else "NOT ACCEPTED")
if (!all(tab$pass)) { print(tab[pass == FALSE]); quit(status = 1L) }
