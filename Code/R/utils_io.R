# =============================================================================
# utils_io.R
# Paths, run tags, config hashing, qs2 caching, seed bookkeeping, logging.
#
# Conventions (Section 5 replication package):
#   * every cached object embeds its config, seeds, package versions, timestamp;
#   * run tags are generated programmatically from the config (never typed by
#     hand), with "." -> "p" so alpha = 0.5 tags as "a0p5", never "a05";
#   * all seeds used anywhere are logged to Results/seeds.csv.
# =============================================================================

library(data.table)
library(here)
library(qs2)

# --- Paths -------------------------------------------------------------------

proj_path <- function(...) here::here(...)

p_data <- function(exp, ...) {
  d <- proj_path("Data", exp)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}

p_results <- function(sub, ...) {
  d <- proj_path("Results", sub)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}

#' Figure path: Results/Figures/<run label>/<experiment>/<file>. Each
#' configuration gets its own self-contained tree.
p_figures <- function(label, exp, ...) {
  d <- proj_path("Results", "Figures", label, exp)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}

# --- Run label ------------------------------------------------------------------
# Names the output tree of one configuration: explicit `label=` wins; else the
# profile, extended by compact tokens for any CLI overrides
# (e.g. "pilot_Kstar20_W5000"); else just the profile.

.LABEL_KEY <- c(K_true = "Kstar", K_grid = "k", J_train = "J", J_eval = "Jev",
                W = "W", alpha_DGP = "a", beta_DGP = "b",
                alpha = "fa", beta = "fb", L = "L", c = "c",
                n_starts = "st", S = "S", S_train = "St")

run_label <- function(cfg, overrides = list(), explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit)) return(explicit)
  if (!length(overrides)) return(cfg$profile)
  toks <- vapply(names(overrides), function(k) {
    v <- gsub("[:,]", "-", overrides[[k]])
    paste0(if (k %in% names(.LABEL_KEY)) .LABEL_KEY[[k]] else k, num2tag(v))
  }, "")
  paste(c(cfg$profile, toks), collapse = "_")
}

# --- Run tags and config hashes ---------------------------------------------

num2tag <- function(x) gsub("\\.", "p", as.character(x))

#' Faithful, filename-safe encoding of a K grid. Distinguishes grids with the
#' same endpoints (10:100 vs c(10,20,40,...)) so result files never collide:
#'   contiguous     -> "k2-20"           (legacy form, existing tags stay valid)
#'   regular step   -> "k5-100by5"
#'   irregular <= 6 -> "k10.20.40.50.70.100"
#'   irregular >  6 -> "k<min>-<max>n<count>x<hash4>"
.k_grid_tag <- function(K_grid) {
  K <- sort(unique(as.integer(K_grid)))
  if (length(K) == 1L) return(paste0("k", K))
  d <- unique(diff(K))
  if (length(d) == 1L) {
    return(if (d == 1L) paste0("k", K[1L], "-", K[length(K)])
           else paste0("k", K[1L], "-", K[length(K)], "by", d))
  }
  if (length(K) <= 6L) return(paste0("k", paste(K, collapse = ".")))
  paste0("k", K[1L], "-", K[length(K)], "n", length(K), "x",
         cfg_hash(K, 4L))
}

run_tag <- function(cfg) {
  stopifnot(!is.null(cfg$experiment), !is.null(cfg$profile))
  parts <- c(
    cfg$experiment,
    cfg$profile,
    if (!is.null(cfg$K_true))  paste0("Kstar", cfg$K_true),
    if (!is.null(cfg$J_train)) paste0("J", cfg$J_train),
    if (!is.null(cfg$W))       paste0("W", cfg$W),
    if (!is.null(cfg$alpha_DGP)) paste0("a", num2tag(cfg$alpha_DGP)),
    if (!is.null(cfg$beta_DGP))  paste0("b", num2tag(cfg$beta_DGP)),
    if (!is.null(cfg$K_grid))  .k_grid_tag(cfg$K_grid),
    # non-default estimator gets its own token so e.g. WarpLDA runs never
    # overwrite VEM result files
    if (!is.null(cfg$fit_method) && cfg$fit_method != "VEM")
      tolower(cfg$fit_method),
    # model/fit priors: append a token ONLY when set, so runs that leave them at
    # the engine default keep their legacy tag/filenames unchanged
    if (!is.null(cfg$alpha) && !is.na(cfg$alpha)) paste0("fa", num2tag(cfg$alpha)),
    if (!is.null(cfg$beta)  && !is.na(cfg$beta))  paste0("fb", num2tag(cfg$beta))
  )
  paste(parts, collapse = "_")
}

cfg_hash <- function(cfg, n = 8L) {
  substr(digest::digest(cfg, algo = "xxhash64"), 1L, n)
}

cache_path <- function(cfg, name) {
  p_data(cfg$experiment, sprintf("%s_%s_%s.qs2", name, run_tag(cfg), cfg_hash(cfg)))
}

# --- qs2 cache with embedded metadata ----------------------------------------

run_meta <- function(cfg) {
  list(
    config     = cfg,
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    r_version  = R.version.string,
    packages   = sapply(
      c("OpTop", "NLPstudio", "topicmodels", "quanteda", "data.table"),
      function(p) as.character(utils::packageVersion(p))
    )
  )
}

cache_put <- function(obj, path, cfg = NULL) {
  attr(obj, "run_meta") <- run_meta(cfg)
  qs2::qs_save(obj, path)
  invisible(path)
}

cache_get <- function(path) qs2::qs_read(path)

#' Compute `expr` unless `path` already exists (incremental re-runs).
with_cache <- function(path, cfg, expr, force = FALSE) {
  if (!force && file.exists(path)) {
    log_msg("cache hit: %s", basename(path))
    return(cache_get(path))
  }
  obj <- force(expr)
  cache_put(obj, path, cfg)
  log_msg("cache write: %s", basename(path))
  obj
}

# --- Seeds -------------------------------------------------------------------
# Deterministic seed schedule per (experiment, scenario, replicate):
# large co-prime strides keep streams disjoint across roles.

make_seeds <- function(seed_base, replicate, scenario_id = 0L) {
  o <- seed_base + 1000000L * scenario_id
  list(
    dgp_seed      = o + 101L  * replicate,
    split_seed    = o + 50021L  + 101L * replicate,
    fit_seed_base = o + 100003L + 1013L * replicate,
    eval_seed     = o + 200003L + 7919L * replicate  # base for E2/E4 eval streams
  )
}

log_seeds <- function(experiment, scenario, replicate, seeds) {
  dt <- data.table(
    experiment = experiment, scenario = scenario, replicate = replicate,
    dgp_seed = seeds$dgp_seed, split_seed = seeds$split_seed,
    fit_seed_base = seeds$fit_seed_base, eval_seed = seeds$eval_seed
  )
  f <- p_results("csv", "seeds.csv")
  if (file.exists(f)) {
    old <- fread(f)
    dt  <- unique(rbindlist(list(old, dt), use.names = TRUE))
  }
  fwrite(dt, f)
  invisible(dt)
}

# --- Tidy result writers ------------------------------------------------------

write_result <- function(dt, name, cfg) {
  stopifnot(is.data.table(dt))
  f <- p_results("csv", sprintf("%s_%s.csv", name, run_tag(cfg)))
  fwrite(dt, f)
  invisible(f)
}

# --- Logging ------------------------------------------------------------------

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
  flush.console()
}
