# =============================================================================
# source_all.R
# Load all Section 5 simulation modules and shared options in one place.
# Usage from any driver:  source(here::here("Code", "R", "source_all.R"))
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(Matrix)
  library(quanteda)
  library(topicmodels)
  library(OpTop)
  library(NLPstudio)
  library(ggplot2)
  library(patchwork)
  library(future.apply)
  library(qs2)
})

# Evaluation-only runs (the revision batch) export OPTOP_NO_FIT=1: a fit-cache
# miss then stops the run instead of estimating a model (see prefit_pool()).
if (nzchar(Sys.getenv("OPTOP_NO_FIT"))) options(optop.no_fit = TRUE)

# explicit seeds are set inside every stochastic helper; silence future's
# RNG-misuse heuristic (fits use topicmodels' own seeded RNG)
options(future.rng.onMisuse = "ignore",
        future.globals.maxSize = 2 * 1024^3)

for (f in c("utils_io.R", "theme_paper.R", "utils_dgp.R", "utils_fit.R",
            "utils_heldout.R", "utils_inference.R", "utils_moment_tests.R",
            "utils_comparators.R", "utils_revision.R")) {
  source(here::here("Code", "R", f))
}

#' Uniform CLI parsing for every entry point.
#'   * `key=value` arguments are configuration overrides (see
#'     apply_overrides() in config/configs.R); `label=name` names the run;
#'   * a bare integer is the worker count;
#'   * remaining words are, in order, the profile and (run_E1 only) the
#'     experiment variant.
#' Examples:
#'   Rscript Code/run_E1.R pilot 8
#'   Rscript Code/run_E1.R full E1b 8 K_true=25 K_grid=5:50
#'   Rscript Code/make_all.R full 10 W=5000 L=500 label=smallvocab
parse_cli <- function(args, default_profile = Sys.getenv("OPTOP_PROFILE",
                                                         "pilot")) {
  is_kv <- grepl("=", args, fixed = TRUE)
  is_num <- !is_kv & grepl("^[0-9]+$", args)
  words <- args[!is_kv & !is_num]

  overrides <- list()
  for (a in args[is_kv]) {
    overrides[[sub("=.*$", "", a)]] <- sub("^[^=]*=", "", a)
  }
  label <- overrides[["label"]]
  overrides[["label"]] <- NULL

  list(
    profile = if (length(words) >= 1L) words[[1L]] else default_profile,
    experiment = if (length(words) >= 2L) words[[2L]] else NULL,
    workers = if (any(is_num)) as.integer(args[is_num][[1L]]) else NULL,
    overrides = overrides,
    label = label
  )
}

#' Resolve a driver's configuration.
#'   cfg_from=<results.qs2>  take the config STORED in an earlier result object
#'       instead of rebuilding it from CLI overrides. Re-scoring runs use this so
#'       the corpus signature, seeds and fit spec -- hence every fit-cache key --
#'       are exactly those of the run that produced the cached models;
#'   out_suffix=<_tag>       appended to every output file name, so a re-run
#'       never overwrites the object it was configured from.
#' Both keys are stripped before the remaining overrides are applied.
resolve_cfg <- function(cli, experiment, profile) {
  ov <- cli$overrides
  cfg_from <- ov$cfg_from; out_suffix <- ov$out_suffix
  ov$cfg_from <- NULL; ov$out_suffix <- NULL
  cfg <- if (!is.null(cfg_from)) {
    if (!file.exists(cfg_from)) stop("cfg_from not found: ", cfg_from, call. = FALSE)
    c0 <- qs2::qs_read(cfg_from)$config
    if (is.null(c0)) stop("no $config in ", cfg_from, call. = FALSE)
    log_msg("config taken from %s", basename(cfg_from))
    c0
  } else get_config(experiment, profile)
  options(optop.revision_suffix = out_suffix %||% "")
  list(cfg = apply_overrides(cfg, ov, strict = FALSE), overrides = ov,
       out_suffix = if (is.null(out_suffix)) "" else out_suffix,
       from_cache = !is.null(cfg_from))
}

setup_parallel <- function(workers = NULL) {
  available <- suppressWarnings(parallel::detectCores(logical = TRUE))
  # detectCores() can be NA in restricted/containerized R sessions.  A serial
  # fallback is safer than allowing min()/if() below to propagate NA and abort
  # even an explicitly requested one-worker smoke run.
  if (length(available) != 1L || !is.finite(available) || available < 1L)
    available <- 1L
  available <- as.integer(available)
  if (is.null(workers)) {
    requested <- suppressWarnings(as.integer(Sys.getenv("OPTOP_WORKERS", NA)))
    workers <- if (is.finite(requested) && requested >= 1L) requested else
      max(1L, available - 1L)
  }
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L)
    stop("workers must be one positive integer", call. = FALSE)
  workers <- min(as.integer(workers), max(1L, available - 1L))
  if (workers > 1L && .Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = workers)
  } else if (workers > 1L) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  log_msg("parallel plan: %s with %d workers",
          class(future::plan())[1], workers)
  workers
}

# --- Shared fit cache: ONE FILE PER (corpus, K) -----------------------------------
# Per-K granularity lets (i) all (seed x K) jobs across replicates pool into a
# single parallel pass with no per-seed barrier, and (ii) pilot fits be reused
# verbatim by the full profile (fit seeds depend on K, not on the grid).

fit_cache_path_k <- function(dgp_signature, K, method, n_starts,
                             fit_seed_base) {
  d <- proj_path("Data", "FITS")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  spec <- list(dgp = dgp_signature, K = K, method = method,
               n_starts = n_starts, fit_seed_base = fit_seed_base)
  file.path(d, sprintf("fit_%s_K%d.qs2", cfg_hash(spec, 12L), K))
}

#' Worker-side real-time scoring log. future buffers a worker's console output
#' until its future resolves, so log_msg() from inside future_lapply appears
#' only at replicate boundaries; file appends escape that (same idiom as
#' prefit.log). Monitor with: tail -f Data/scoring.log
scoring_log <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)),
      file = proj_path("Data", "scoring.log"), append = TRUE)
}

#' Worker for prefit_pool(): TOP-LEVEL on purpose -- a closure defined inside
#' prefit_pool() would drag the whole job list into the exported environment.
#' Appends one line per completed job to Data/FITS/prefit.log (tail -f it to
#' monitor long runs) and ticks the progressr progressor when one is passed.
.prefit_worker <- function(jb, p = NULL) {
  t0 <- proc.time()[["elapsed"]]
  out <- .fit_one_k(jb$dtm, jb$K, jb$method, jb$n_starts, jb$seed_base,
                    fit_alpha = jb$fit_alpha, fit_beta = jb$fit_beta)
  qs2::qs_save(out, jb$path)
  cat(sprintf("[%s] done K=%d %s (%d start%s) in %.1fs\n",
              format(Sys.time(), "%H:%M:%S"), jb$K, jb$method, jb$n_starts,
              if (jb$n_starts > 1L) "s" else "",
              proc.time()[["elapsed"]] - t0),
      file = file.path(dirname(jb$path), "prefit.log"), append = TRUE)
  if (!is.null(p)) p(sprintf("K=%d %s", jb$K, jb$method))
  jb$path
}

#' Pooled prefit: fit every missing (corpus, K) job of a batch in ONE parallel
#' pass. `jobs` is a list of lists with fields: dtm, K_grid, method, n_starts,
#' fit_seed_base, signature.
prefit_pool <- function(jobs) {
  todo <- list()
  for (jb in jobs) {
    for (K in jb$K_grid) {
      p <- fit_cache_path_k(jb$signature, K, jb$method, jb$n_starts,
                            jb$fit_seed_base)
      if (!file.exists(p)) {
        fit_alpha <- if (!is.null(jb$signature$fit_alpha))
          jb$signature$fit_alpha else NA_real_
        fit_beta <- if (!is.null(jb$signature$fit_beta))
          jb$signature$fit_beta else NA_real_
        todo[[length(todo) + 1L]] <- list(
          dtm = jb$dtm, K = K, method = jb$method,
          n_starts = jb$n_starts, seed_base = jb$fit_seed_base, path = p,
          fit_alpha = fit_alpha, fit_beta = fit_beta)
      }
    }
  }
  if (!length(todo)) { log_msg("fit pool: all %d grids cached", length(jobs)); return(invisible(0L)) }
  # No-fit guard: evaluation-only scripts (postprocessing, rescoring, the
  # revision batch) set options(optop.no_fit = TRUE). A cache miss then means
  # the corpus signature or fit spec drifted from the run that produced the
  # cached models, and must never be "repaired" by silently estimating new ones.
  if (isTRUE(getOption("optop.no_fit", FALSE))) {
    miss <- vapply(todo, function(jb) basename(jb$path), "")
    stop(sprintf(paste0("optop.no_fit is set but %d fit(s) are missing from ",
                        "Data/FITS (first: %s). Refusing to estimate models; ",
                        "check the corpus signature / fit spec."),
                 length(miss), miss[1L]), call. = FALSE)
  }
  # longest jobs first so the parallel tail stays short
  todo <- todo[order(-vapply(todo, `[[`, 1L, "K"))]

  # upfront composition so slow arms are visible BEFORE committing hours
  bd <- data.table(
    K = vapply(todo, `[[`, 1L, "K"),
    method = vapply(todo, `[[`, "", "method"),
    n_starts = vapply(todo, `[[`, 1L, "n_starts")
  )
  bs <- bd[, .(n = .N, kmin = min(K), kmax = max(K), starts = max(n_starts)),
           by = method]
  log_msg("fit pool: %d (corpus x K) jobs, longest K first -- %s",
          length(todo),
          paste(sprintf("%s %d (K %d-%d, %d start%s)", bs$method, bs$n,
                        bs$kmin, bs$kmax, bs$starts,
                        ifelse(bs$starts > 1L, "s", "")), collapse = ", "))
  log_msg("progress: tail -f %s",
          proj_path("Data", "FITS", "prefit.log"))
  big_J <- max(vapply(jobs, function(jb) nrow(jb$dtm), 0L))
  slow <- bd$method %in% c("Gibbs", "VEM")
  if (any(slow) && big_J >= 2000L) {
    log_msg(paste0("NOTE: %d %s job(s) on corpora with up to J=%d documents ",
                   "-- these can take many minutes to ~an hour EACH ",
                   "(fit_method=WarpLDA is ~100x faster)"),
            sum(slow), paste(unique(bd$method[slow]), collapse = "/"), big_J)
  }

  # one future per job (future.scheduling = Inf) so the bar ticks per fit;
  # the file log is the dependable channel if the handler cannot render
  res <- tryCatch(
    progressr::with_progress({
      p <- progressr::progressor(steps = length(todo))
      future_lapply(todo, .prefit_worker, p = p,
                    future.seed = NULL, future.scheduling = Inf)
    }),
    error = function(e) future_lapply(todo, .prefit_worker,
                                      future.seed = NULL,
                                      future.scheduling = Inf)
  )
  invisible(res)
}

#' Fetch a cached LDA grid, fitting any missing K (serially-parallel fallback).
get_fits_cached <- function(dtm, K_grid, method, n_starts, fit_seed_base,
                            dgp_signature, force = FALSE) {
  paths <- vapply(K_grid, function(K)
    fit_cache_path_k(dgp_signature, K, method, n_starts, fit_seed_base), "")
  missing <- K_grid[!file.exists(paths)]
  if (force) missing <- K_grid
  if (length(missing)) {
    prefit_pool(list(list(dtm = dtm, K_grid = missing, method = method,
                          n_starts = n_starts, fit_seed_base = fit_seed_base,
                          signature = dgp_signature)))
  }
  outs <- lapply(paths, qs2::qs_read)
  models <- lapply(outs, `[[`, "fit")
  names(models) <- as.character(K_grid)
  list(models = models,
       diagnostics = rbindlist(lapply(outs, `[[`, "diag"), fill = TRUE))
}
