# =============================================================================
# list_runs.R -- discover cached experiment results and how to re-plot them.
#
#   Rscript Code/list_runs.R [filter ...]
#
# Scans Data/E*/*_results_*.qs2 (the caches make_figures.R / make_tables.R read)
# and prints, per cached run, a ready-to-paste make_figures.R command.
#
# Why this exists: data files are named by a parameter TAG (run_tag(): K_true,
# J_train, W, alpha_DGP, beta_DGP, K_grid, fit_method, optional fit priors) --
# never by the `label=` you gave
# make_all.R, and never by S / J_eval / gibbs_S. So to re-plot a finished run you
# must reproduce its tag fields. This tool prints exactly that command, so you
# never have to decode a filename by hand.
#
# The command's tag fields are read from the FILENAME (ground truth for what is
# on disk); the embedded config supplies only the non-tag metadata (S, J_eval,
# n_starts). Each command is self-checked to reconstruct its own cache file.
#
# Optional args are case-insensitive substrings ANDed against the tag:
#   Rscript Code/list_runs.R E1            # only E1 runs
#   Rscript Code/list_runs.R full warplda  # full-profile WarpLDA runs
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
suppressWarnings(suppressMessages({
  source(here::here("Code", "R", "source_all.R"))
  source(here::here("Code", "config", "configs.R"))
}))

filters <- tolower(args)

# --- helpers ----------------------------------------------------------------

# tag-bearing config fields, in run_tag() order (see utils_io.R:82).  `alpha`
# and `beta` are optional fitted-model priors; the DGP priors have explicit
# names so a rerun cannot silently swap simulation and fitting semantics.
.TAG_FIELDS   <- c("K_true", "J_train", "W", "alpha_DGP", "beta_DGP",
                   "K_grid", "fit_method", "alpha", "beta")
.METHOD_CANON <- c(warplda = "WarpLDA", gibbs = "Gibbs", vem = "VEM")

# inverse of the CLI grid syntax accepted by .parse_override_value():
#   length 1 -> "8"   regular step d -> "a:b" (d==1) or "a:b:d"   else "a,b,c"
grid_to_cli <- function(K) {
  K <- sort(unique(as.integer(K)))
  if (!length(K)) return("")
  if (length(K) == 1L) return(as.character(K))
  d <- unique(diff(K))
  if (length(d) == 1L)
    return(if (d == 1L) sprintf("%d:%d", K[1L], K[length(K)])
           else          sprintf("%d:%d:%d", K[1L], K[length(K)], d))
  paste(K, collapse = ",")
}

# a K-grid tag token (k2-20 / k5-100by5 / k10.20.40 / k8) back to an integer
# vector; the hashed irregular form (k<lo>-<hi>n<N>x<hash>) is unrecoverable -> NULL
k_tag_to_grid <- function(tok) {
  s <- sub("^k", "", tok)
  if (grepl("x", s, fixed = TRUE)) return(NULL)
  if (grepl(".", s, fixed = TRUE))
    return(as.integer(strsplit(s, ".", fixed = TRUE)[[1]]))
  if (grepl("by", s, fixed = TRUE)) {
    ab  <- strsplit(s, "by", fixed = TRUE)[[1]]
    rng <- as.integer(strsplit(ab[1L], "-", fixed = TRUE)[[1]])
    return(seq.int(rng[1L], rng[2L], as.integer(ab[2L])))
  }
  if (grepl("-", s, fixed = TRUE)) {
    rng <- as.integer(strsplit(s, "-", fixed = TRUE)[[1]])
    return(seq.int(rng[1L], rng[2L]))
  }
  as.integer(s)
}

# rebuild the tag fields from a filename tag -- authoritative for what run_tag
# produced, hence for which file exists on disk (the embedded config can carry a
# different K_grid, e.g. E5's grid-extension arm).
parse_tag <- function(tag) {
  p   <- strsplit(tag, "_", fixed = TRUE)[[1]]
  cfg <- list(experiment = p[1L], profile = p[2L])
  for (tok in p[-(1:2)]) {
    if      (grepl("^Kstar", tok))   cfg$K_true  <- as.integer(sub("^Kstar", "", tok))
    else if (grepl("^J[0-9]", tok))  cfg$J_train <- as.integer(sub("^J", "", tok))
    else if (grepl("^W[0-9]", tok))  cfg$W       <- as.integer(sub("^W", "", tok))
    else if (grepl("^fa[0-9p]", tok)) cfg$alpha <- as.numeric(gsub("p", ".", sub("^fa", "", tok)))
    else if (grepl("^fb[0-9p]", tok)) cfg$beta  <- as.numeric(gsub("p", ".", sub("^fb", "", tok)))
    else if (grepl("^a[0-9p]", tok)) cfg$alpha_DGP <- as.numeric(gsub("p", ".", sub("^a", "", tok)))
    else if (grepl("^b[0-9p]", tok)) cfg$beta_DGP  <- as.numeric(gsub("p", ".", sub("^b", "", tok)))
    else if (grepl("^k[0-9]", tok))  cfg$K_grid  <- k_tag_to_grid(tok)
    else { m <- .METHOD_CANON[tolower(tok)]; cfg$fit_method <- if (is.na(m)) tok else unname(m) }
  }
  cfg
}

# the run_tag() token a single field would contribute -- lets us diff a field
# against the profile default regardless of numeric type (10L vs 10, 2:20 vs
# 2:20L), so the override set stays minimal and never spuriously verbose.
tok_of <- function(f, v) {
  if (is.null(v)) return(NA_character_)
  switch(f,
    K_true     = paste0("Kstar", v),
    J_train    = paste0("J", v),
    W          = paste0("W", v),
    alpha_DGP  = paste0("a", num2tag(v)),
    beta_DGP   = paste0("b", num2tag(v)),
    K_grid     = .k_grid_tag(v),
    fit_method = if (identical(as.character(v), "VEM")) "" else tolower(as.character(v)),
    alpha      = paste0("fa", num2tag(v)),
    beta       = paste0("fb", num2tag(v)),
    as.character(v))
}

# minimal CLI overrides: tag fields whose tag token differs from the profile
# default. Named list of string values (as apply_overrides expects).
overrides_for <- function(cfg) {
  base <- tryCatch(get_config(cfg$experiment, cfg$profile), error = function(e) NULL)
  ov <- list()
  for (f in .TAG_FIELDS) {
    v <- cfg[[f]]
    if (is.null(v)) next
    bv <- if (!is.null(base)) base[[f]] else NULL
    if (!is.null(bv) && identical(tok_of(f, v), tok_of(f, bv))) next
    ov[[f]] <- if (f == "K_grid") grid_to_cli(v) else as.character(v)
  }
  ov
}

# tag encoded in a cache filename: strip "<exp>_results_" prefix and ".qs2"
tag_of <- function(path)
  sub("\\.qs2$", "", sub("^[a-z0-9]+_results_", "", basename(path)))

# embedded config (for non-tag metadata only): out$config, else the run_meta attr
read_cfg <- function(path) {
  obj <- tryCatch(cache_get(path), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  cfg <- obj$config
  if (is.null(cfg)) cfg <- attr(obj, "run_meta")$config
  cfg
}

human_size <- function(b) {
  u <- c("B", "KB", "MB", "GB"); i <- 1L
  while (b >= 1024 && i < length(u)) { b <- b / 1024; i <- i + 1L }
  sprintf("%.1f %s", b, u[i])
}

# --- scan -------------------------------------------------------------------

files <- Sys.glob(file.path(proj_path("Data"), "E*", "*_results_*.qs2"))
if (length(filters))
  files <- files[vapply(files, function(f) {
    t <- tolower(tag_of(f))
    all(vapply(filters, function(q) grepl(q, t, fixed = TRUE), logical(1)))
  }, logical(1))]

if (!length(files)) {
  cat(if (length(args))
        sprintf("No cached runs match: %s\n", paste(args, collapse = " "))
      else "No cached runs found under Data/E*/.\n")
  quit(save = "no", status = 0L)
}

# group by experiment (filename prefix e1_/e1b_/e2_...), newest first within each
mtime <- file.info(files)$mtime
files <- files[order(basename(files), -as.numeric(mtime))]

cat("Cached experiment runs  --  Data/E*/*_results_*.qs2\n")
cat("Each command regenerates that run's figures from its cache (no recompute).\n")
cat("  - data is keyed by parameters (the tag); label= only names the OUTPUT folder\n")
cat("  - append  label=NAME  to pick Results/Figures/NAME/  (else auto-named from params)\n")
cat("  - swap  make_figures -> make_tables  for that run's tables + Excel sheet\n")
cat(strrep("-", 90), "\n", sep = "")

n <- 0L
for (f in files) {
  tag  <- tag_of(f)
  tcfg <- parse_tag(tag)                 # tag fields, authoritative (from filename)
  mcfg <- read_cfg(f)                    # embedded config, metadata only
  if (is.null(tcfg$K_grid) && !is.null(mcfg))   # hashed grid unrecoverable from name
    tcfg$K_grid <- mcfg$K_grid
  fi <- file.info(f)

  hdr <- sprintf("[%s / %s]  K*=%s  J_train=%s  W=%s  grid=%s  %s",
                 tcfg$experiment, tcfg$profile,
                 tcfg$K_true %||% "?", tcfg$J_train %||% "?", tcfg$W %||% "?",
                 if (!is.null(tcfg$K_grid)) grid_to_cli(tcfg$K_grid) else "?",
                 tcfg$fit_method %||% "VEM")
  extra <- character(0)
  if (!is.null(tcfg$alpha_DGP) && !isTRUE(all.equal(tcfg$alpha_DGP, 0.5)))
    extra <- c(extra, paste0("alpha_DGP=", tcfg$alpha_DGP))
  if (!is.null(tcfg$beta_DGP) && !isTRUE(all.equal(tcfg$beta_DGP, 0.01)))
    extra <- c(extra, paste0("beta_DGP=", tcfg$beta_DGP))
  # `$.` partially matches by default in R: without [[ ]] an absent fitted
  # `beta` would be mistaken for `beta_DGP` and falsely displayed as an
  # explicit fitted-model prior.
  fit_alpha <- tcfg[["alpha"]]
  fit_beta  <- tcfg[["beta"]]
  if (!is.null(fit_alpha)) extra <- c(extra, paste0("alpha=", fit_alpha))
  if (!is.null(fit_beta))  extra <- c(extra, paste0("beta=", fit_beta))
  if (length(extra)) hdr <- paste0(hdr, "  ", paste(extra, collapse = " "))
  cat("\n", hdr, "\n", sep = "")

  meta <- c(human_size(fi$size), format(fi$mtime, "%Y-%m-%d %H:%M"))
  if (!is.null(mcfg)) {                  # S / J_eval / n_starts live only in the config
    for (k in c("S", "S_train", "n_starts", "J_eval", "gibbs_S"))
      if (!is.null(mcfg[[k]])) meta <- c(meta, sprintf("%s=%s", k, mcfg[[k]]))
  } else meta <- c(meta, "(no embedded config)")
  cat("    ", paste(meta, collapse = "  |  "), "\n", sep = "")

  ov  <- overrides_for(tcfg)
  cmd <- sprintf("Rscript Code/make_figures.R %s exp=%s%s",
                 tcfg$profile, tcfg$experiment,
                 if (length(ov)) paste0(" ", paste(sprintf("%s=%s", names(ov), unlist(ov)),
                                                    collapse = " ")) else "")
  # self-check: the emitted overrides must rebuild this exact cache tag
  recon <- tryCatch(run_tag(apply_overrides(get_config(tcfg$experiment, tcfg$profile),
                                            ov, strict = FALSE)),
                    error = function(e) NA_character_)
  if (!identical(recon, tag))
    cmd <- paste0(cmd, sprintf("   # NOTE: check K_grid -- reconstructs '%s'", recon %||% "NA"))
  cat("    ", cmd, "\n", sep = "")
  n <- n + 1L
}
cat(strrep("-", 90), "\n", sep = "")
cat(sprintf("%d cached run%s. Re-run with a filter (e.g. `E1`, `full`, `warplda`) to narrow.\n",
            n, if (n == 1L) "" else "s"))
