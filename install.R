# ------------------------------------------------------------------------------
# install.R -- installs the package versions used for the published results.
# Run once from the package root:  Rscript install.R
#
# CRAN packages are installed only if absent or older than the pinned version
# (newer versions are accepted; the exact environments are in sessionInfo.txt).
# OpTop and NLPstudio are pinned to the exact GitHub commits of the environment
# that produced the September 2026 robustness arms and regenerated every
# shipped exhibit: OpTop 0.20.0 and NLPstudio 1.2.0. The July 2026 baseline
# caches were produced under OpTop 0.14.0/0.14.1 + NLPstudio 1.1.1; the
# pipeline runs unchanged under both (README.md, Section 2).
# ------------------------------------------------------------------------------

if (getRversion() < "4.6") stop("R >= 4.6 is required (shipped results: R 4.6.1).")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran_pins <- c(
  data.table   = "1.18.4",
  Matrix       = "1.7.5",
  future.apply = "1.20.2",
  ggplot2      = "4.0.3",   # >= 4.0 required (geom_errorbar orientation API)
  patchwork    = "1.3.2",
  qs2          = "0.2.2",
  tinytable    = "0.17.0",
  writexl      = "1.5.4",
  here         = "1.0.2",
  digest       = "0.6.39",
  MASS         = "7.3.65",
  topicmodels  = "0.2.17",
  quanteda     = "4.4",
  text2vec     = "0.6.6"    # provides the WarpLDA engine
)

for (p in names(cran_pins)) {
  ok <- requireNamespace(p, quietly = TRUE) &&
    packageVersion(p) >= package_version(cran_pins[[p]])
  if (!ok) {
    message("installing ", p, " (>= ", cran_pins[[p]], ")")
    install.packages(p)
  }
}

# --- GitHub pins (exact commits behind the published results) ----------------
if (!requireNamespace("OpTop", quietly = TRUE) ||
    packageVersion("OpTop") < "0.20.0") {
  remotes::install_github("contefranz/OpTop@1166daed33")       # v0.20.0 (main, 2026-07-22)
}
if (!requireNamespace("NLPstudio", quietly = TRUE) ||
    packageVersion("NLPstudio") < "1.2.0") {
  remotes::install_github("contefranz/NLPstudio@771b42e328")   # v1.2.0 (main, 2026-07-23)
}
message("OpTop ", packageVersion("OpTop"), " | NLPstudio ", packageVersion("NLPstudio"),
        " | R ", getRversion())

message("\nAll dependencies satisfied. Next steps:")
message("  Rscript Code/tests_unit.R     # correctness gates")
message("  ./reproduce.sh exhibits       # rebuild every exhibit from the caches")
