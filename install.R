# ------------------------------------------------------------------------------
# install.R -- installs the exact software environment behind every shipped result.
#
#   Rscript install.R              # install into the default R library
#   Rscript install.R lib=.Rlib    # recommended: a project-local library (.Rlib/,
#                                  # git-ignored); reproduce.sh uses it automatically
#
# CRAN packages come from the Posit Package Manager snapshot of 2026-07-01, which
# serves exactly the versions that produced the shipped results (R 4.6.1; all 70
# packages of the dependency closure were checked against it on 9 Oct 2026; see
# sessionInfo.txt). OpTop 0.20.1 and NLPstudio 1.2.0 are installed from their pinned
# GitHub commits; an installed copy is accepted only if its RemoteSha matches.
#
# OpTop is compiled from source (C++ with RcppArmadillo and OpenMP): it needs a C++17
# toolchain and a Fortran runtime (macOS: Xcode Command Line Tools + the CRAN gfortran
# build; Linux: build-essential and gfortran; Windows: Rtools). Snapshot packages are
# installed as binaries where the snapshot offers them for your platform, else from
# source. Unit gate U20 needs OpTop 0.20.1 exactly (its classed word-null warning).
# ------------------------------------------------------------------------------

SNAPSHOT <- "https://packagemanager.posit.co/cran/2026-07-01"
options(repos = c(CRAN = SNAPSHOT))
nc <- parallel::detectCores()                       # parallel package builds; changes speed only
options(Ncpus = if (is.na(nc)) 1L else max(1L, nc - 1L))
if (getRversion() < "4.6") stop("R >= 4.6 is required (shipped results: R 4.6.1).")

args <- commandArgs(trailingOnly = TRUE)
lib_arg <- sub("^lib=", "", grep("^lib=", args, value = TRUE))
if (length(lib_arg)) {
  dir.create(lib_arg, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(normalizePath(lib_arg), .libPaths()))
  message("installing into the project library ", normalizePath(lib_arg))
}
lib <- .libPaths()[1L]

# --- CRAN packages used directly by the replication code (exact versions) -------
pins <- c(
  data.table   = "1.18.4",  Matrix      = "1.7-5",   future       = "1.70.0",
  future.apply = "1.20.2",  ggplot2     = "4.0.3",   patchwork    = "1.3.2",
  qs2          = "0.2.2",   tinytable   = "0.17.0",  writexl      = "1.5.4",
  here         = "1.0.2",   digest      = "0.6.39",  MASS         = "7.3-65",
  topicmodels  = "0.2-17",  quanteda    = "4.4",     text2vec     = "0.6.6",   # WarpLDA engine
  jsonlite     = "2.0.0",   remotes     = "2.5.0",   progressr    = "0.19.0",
  RhpcBLASctl  = "0.23-42", ragg       = "1.5.2"    # ragg: the PNG device ggsave uses (pixel-identical PNGs)
)
has <- function(p, v) {
  ok <- suppressWarnings(requireNamespace(p, quietly = TRUE))
  ok && packageVersion(p) == package_version(v)
}
need <- names(pins)[!mapply(has, names(pins), pins)]
if (length(need)) {
  message("installing from the 2026-07-01 snapshot: ", paste(need, collapse = ", "))
  install.packages(need, lib = lib)
}
bad <- names(pins)[!mapply(has, names(pins), pins)]
if (length(bad)) stop("version mismatch after installation: ",
  paste(sprintf("%s %s (expected %s)", bad,
    vapply(bad, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "missing"), ""),
    pins[bad]), collapse = "; "),
  "\nA package loaded by another R process can block an update: close R sessions and rerun,",
  "\nor install into a fresh library with  Rscript install.R lib=.Rlib", call. = FALSE)

# --- GitHub packages, pinned to exact commits ------------------------------------
gh <- c(OpTop     = "contefranz/OpTop@cba1273ea8061d3afb7cb26344779d95c492c90f",       # v0.20.1 (2026-09-19)
        NLPstudio = "contefranz/NLPstudio@771b42e328098fae6e5d9ed98e779b130793aff9")   # v1.2.0  (2026-07-23)
sha_ok <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE)) &&
  identical(packageDescription(p)$RemoteSha, sub(".*@", "", gh[[p]]))
for (p in names(gh)) if (!sha_ok(p))
  remotes::install_github(gh[[p]], lib = lib, upgrade = "never", force = TRUE, dependencies = NA)
if (!all(vapply(names(gh), sha_ok, NA))) stop("OpTop/NLPstudio are not at the pinned commits.", call. = FALSE)

# --- whole dependency closure against the snapshot (informative) -----------------
ip <- installed.packages()
deps <- unique(c(names(pins), unlist(tools::package_dependencies(c(names(pins), names(gh)), db = ip,
  recursive = TRUE, which = c("Depends", "Imports", "LinkingTo")))))
deps <- setdiff(deps, c(rownames(ip)[ip[, "Priority"] %in% "base"], names(gh)))
ap <- available.packages(type = "source")
inst <- ip[match(deps, ip[, "Package"]), "Version"]; snap <- ap[match(deps, ap[, "Package"]), "Version"]
off <- deps[is.na(snap) | inst != snap]
if (length(off)) {
  message("note: ", length(off), " of ", length(deps), " dependencies differ from the snapshot: ",
          paste(sprintf("%s %s (snapshot %s)", off, inst[match(off, deps)], snap[match(off, deps)]), collapse = ", "),
          "\n      for an exact environment use a fresh library:  Rscript install.R lib=.Rlib")
} else message("all ", length(deps), " dependencies match the 2026-07-01 snapshot")

message("OpTop ", packageVersion("OpTop"), " @ ", substr(packageDescription("OpTop")$RemoteSha, 1, 10),
        " | NLPstudio ", packageVersion("NLPstudio"), " @ ", substr(packageDescription("NLPstudio")$RemoteSha, 1, 10),
        " | R ", getRversion(), " | library ", lib)
message("\nNext steps:")
message("  ./reproduce.sh test             # unit gates U1-U25")
message("  ./reproduce.sh fetch results    # result objects from the Zenodo data record")
message("  ./reproduce.sh paper            # every exhibit of the paper and supplement")
