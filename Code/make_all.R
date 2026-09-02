# =============================================================================
# make_all.R -- orchestrate the full Section 5 pipeline.
# Usage (from the project root):
#   Rscript Code/make_all.R smoke   # end-to-end exercise, minutes
#   Rscript Code/make_all.R pilot   # ~1 h on ~11 workers
#   Rscript Code/make_all.R full    # paper scale (overnight)
# Each stage runs in a fresh R session; fit caches make re-runs incremental.
# =============================================================================

# Usage: Rscript Code/make_all.R [smoke|pilot|full] [workers] [key=value ...] [label=name]
# key=value overrides and label= are forwarded to every stage; each stage applies
# the keys its experiment uses and ignores the rest with a notice (see
# apply_overrides(strict = FALSE) in config/configs.R). A key valid for no
# experiment (a typo) fails fast below before any stage runs.
args <- commandArgs(trailingOnly = TRUE)
is_kv <- grepl("=", args, fixed = TRUE)
is_num <- !is_kv & grepl("^[0-9]+$", args)
words <- args[!is_kv & !is_num]
PROFILE <- if (length(words)) words[[1L]] else "smoke"
stopifnot(PROFILE %in% c("smoke", "pilot", "full"))
if (any(is_num)) {
  Sys.setenv(OPTOP_WORKERS = args[is_num][[1L]])   # inherited by every stage
  cat(sprintf("worker count: %s\n", args[is_num][[1L]]))
}
FORWARD <- args[is_kv]
if (length(FORWARD)) cat("forwarded:", paste(FORWARD, collapse = " "), "\n")

# Fail-fast typo guard: a forwarded key must be valid for at least ONE
# experiment (each stage then applies the keys it uses and ignores the rest via
# apply_overrides(strict = FALSE)). Degrades to no-guard if configs can't load.
if (length(FORWARD)) tryCatch({
  suppressWarnings(suppressMessages(
    source(here::here("Code", "config", "configs.R"))))
  valid <- unique(c("L", "label", "exp", unlist(lapply(
    c("E1", "E1b", "E1c", "E2", "E3", "E4", "E5", "E6"),
    function(e) names(get_config(e, PROFILE))))))
  keys <- sub("=.*$", "", FORWARD)
  bad <- setdiff(keys, valid)
  if (length(bad)) {
    hint <- vapply(bad, function(k)
      valid[which.min(adist(k, valid, ignore.case = TRUE))], "")
    stop(sprintf("unknown override(s): %s\n(did you mean: %s?)\nValid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(sprintf("%s -> %s", bad, hint), collapse = "; "),
                 paste(sort(valid), collapse = ", ")), call. = FALSE)
  }
}, error = function(e) if (grepl("^unknown override", conditionMessage(e)))
  stop(e) else message("override guard skipped: ", conditionMessage(e)))

stage <- function(script, ...) {
  t0 <- Sys.time()
  cat(sprintf("\n========== %s [%s] ==========\n", script, PROFILE))
  status <- system2("Rscript", c(file.path("Code", script), PROFILE, ...,
                                 FORWARD))
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat(sprintf("---------- %s finished in %.1f min (status %d)\n",
              script, dt, status))
  if (status != 0L) stop(script, " failed with status ", status)
  invisible(dt)
}

t_all <- Sys.time()
stage("tests_unit.R")
stage("run_E1.R")
stage("run_E2.R")
stage("run_E3.R")
stage("run_E4.R")
stage("run_E5.R")
stage("run_E6.R")
# E1b/E1c robustness variants are no longer part of the full orchestration
# (~110 core-h, no draft float, no referee item). Run manually if needed:
#   Rscript Code/run_E1.R full E1b <workers>
stage("make_figures.R")
stage("make_tables.R")
cat(sprintf("\n=== ALL DONE [%s] in %.1f min ===\n", PROFILE,
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
