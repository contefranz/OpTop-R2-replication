# =============================================================================
# utils_revision.R -- revision-only provenance, checkpoints and diagnostics.
#
# Nothing here touches a fit-cache key (corpus signature, dgp_signature,
# fit_cache_path_k, run_tag), so every cached topic model stays reachable.
#
#   * revision_provenance()    what defines a SCORE in objects written from now
#                              on: package version / commit and the convention
#                              of the held-out word-level null deviance;
#   * word_null_convention_of() / revision_word_null()
#                              the ONLY route by which a word-level null is
#                              converted between conventions -- legacy inputs
#                              are identified by content hash, new inputs by
#                              their stamp, anything else is an error;
#   * revision_checkpoint()    resumable units of evaluation work, reused only
#                              when configuration, package and scoring code are
#                              identical;
#   * support_resolution()     how coarse the harmonised support is, per
#                              document and candidate.
# =============================================================================

WORD_NULL_CONVENTIONS <- c(
  # OpTop <= 0.20.1 word-level null on an EXTERNAL baseline: 2 sum_j N log(N/B),
  # without the Poisson linear term (0.20.1 warns about it, the kernel is unchanged)
  legacy   = "legacy_log_only",
  # Section 2.4 of the paper, computed by .word_null_dev_poisson() on the scored
  # tokens: 2 sum_j [N log(N/B) - (N - B)], B_jw = L_j^{scored} pi_w^{train}
  poisson  = "poisson_scored_tokens")

revision_provenance <- function() {
  d <- utils::packageDescription("OpTop")
  list(word_null_convention = WORD_NULL_CONVENTIONS[["poisson"]],
       optop_version = as.character(utils::packageVersion("OpTop")),
       optop_sha = d$RemoteSha %||% d$GithubSHA1 %||% NA_character_)
}

# --- which convention does a result object carry? ---------------------------------

#' Pre-revision result objects that contain HELD-OUT word-level deviances written
#' under the legacy convention, identified by the sha256 of the file. A file that
#' is neither stamped (run_meta$scoring) nor listed here is AMBIGUOUS and is
#' refused: the convention is never inferred from a package version, a file name
#' or the sign of an index.
LEGACY_WORD_NULL_INPUTS <- c(
  "e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2" =
    "7c7bc21e03cd6678abaaf0b934fb07c7520c22b332fa4434b5fa962cf04c950d",
  "e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p5_fb0p01.qs2" =
    "9c9f841fe93b8f34d664d08ecaffc7b3e2489a7139ddeb088128590d880b71fe",
  "e6_results_E6_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2" =
    "97c52070fe620a3fab083b277163702f4fa77495418ddfb1d2a8e49910670c33",
  "mdna_results_MDNA_2015_2016.qs2" =
    "ffbeb0be8e828e743cbde539b6a5dc4a4631c7c8b918c2e8e541eedf09bfbb1d",
  "mdna_results_MDNA_2015_2016_c05.qs2" =
    "3886dbb34a0f3ec5f84f585c8fefa98e0e707f3add14e76cee805fe1e9d0e715",
  "mdna_results_MDNA_2015_2016_c2.qs2" =
    "314ea428a5c2a19af4f50095da9a638bc5d5ccfd17fe3953a1c2cd5b3f9be525",
  "mdna_results_MDNA_2015_2016_g100.qs2" =
    "a3ef050b66c80227c2305c29a6bd887756b87a98a142c4da5c7aaade06a2c308")

#' Convention of the held-out word-level null stored in result object `obj`, read
#' from `path`. Stamped objects answer for themselves; unstamped objects must be
#' one of the frozen legacy files above, byte for byte.
word_null_convention_of <- function(obj, path) {
  stamp <- attr(obj, "run_meta")$scoring$word_null_convention %||%
    obj$scoring$word_null_convention
  if (!is.null(stamp)) {
    if (length(stamp) != 1L || !stamp %in% WORD_NULL_CONVENTIONS)
      stop("unknown word-null convention stamp '", paste(stamp, collapse = ","),
           "' in ", path, call. = FALSE)
    return(unname(stamp))
  }
  known <- LEGACY_WORD_NULL_INPUTS[basename(path)]
  if (is.na(known))
    stop("AMBIGUOUS word-null convention: ", path, " carries no scoring stamp ",
         "and is not a registered legacy input. Refusing to guess.", call. = FALSE)
  got <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(got, unname(known)))
    stop("AMBIGUOUS word-null convention: ", path, " has the name of a ",
         "registered legacy input but not its content (sha256 ", got, ").",
         call. = FALSE)
  WORD_NULL_CONVENTIONS[["legacy"]]
}

#' Bring a word table (columns d_model, d_null) to the paper's convention. The
#' conversion is applied ONLY to tables declared legacy; a table already in
#' Poisson form passes through untouched, so the correction can never be
#' applied twice. `N_w`, `B_w`: scored-corpus count and baseline expected count
#' of each word, aligned with the rows of `word`.
revision_word_null <- function(word, N_w, B_w, convention) {
  if (missing(convention) || length(convention) != 1L || is.na(convention) ||
      !convention %in% WORD_NULL_CONVENTIONS)
    stop("Unknown word-null convention", call. = FALSE)
  stopifnot(is.data.table(word), all(c("d_model", "d_null") %in% names(word)),
            nrow(word) == length(N_w), length(N_w) == length(B_w))
  ans <- data.table::copy(word)
  if (convention == WORD_NULL_CONVENTIONS[["legacy"]])
    ans[, d_null := d_null - 2 * (N_w - B_w)]
  ans[, r2_word := ifelse(d_null > 0, 1 - d_model / d_null, NA_real_)]
  ans[]
}

# --- checkpoints ---------------------------------------------------------------------

#' Files whose content determines a SCORE in the running process: the shared
#' modules, the configuration, and the script being run (Rscript --file=). Other
#' drivers and every exhibit / post-processing script are deliberately excluded,
#' so repairing a figure script -- or a different driver -- never invalidates a
#' night of scoring. (revision_batch.py builds each stage's identity the same way.)
revision_code_files <- function() {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  script <- if (length(script) && file.exists(script[1L])) normalizePath(script[1L]) else
    character()
  sort(unique(c(list.files(proj_path("Code", "R"), "\\.R$", full.names = TRUE),
                proj_path("Code", "config", "configs.R"), script)))
}

.revision_memo <- new.env(parent = emptyenv())

#' sha256 over (file name, md5) pairs; names only, so the hash does not depend on
#' where the project lives. Computed once per process (forked workers inherit it).
revision_code_hash <- function(refresh = FALSE) {
  if (refresh || is.null(.revision_memo$code)) {
    f <- revision_code_files()
    md5 <- unname(tools::md5sum(f))
    if (anyNA(md5)) stop("cannot hash scoring code: ", paste(f[is.na(md5)], collapse = ", "))
    .revision_memo$code <- digest::digest(paste(basename(f), md5), algo = "sha256")
  }
  .revision_memo$code
}

revision_identity <- function(cfg) {
  list(config = cfg, provenance = revision_provenance(), code = revision_code_hash())
}

#' Evaluate `expr` once per (suffix, stage, configuration, key) and reuse the
#' stored value afterwards -- but only if configuration, package and scoring code
#' are IDENTICAL to those that produced it. A mismatch is an error, never a
#' silent recomputation over stale neighbours: choose a new out_suffix, or delete
#' the stale checkpoint directory. Without an out_suffix (ordinary runs, unit
#' fixtures) nothing is written and `expr` is simply evaluated.
revision_checkpoint <- function(stage, cfg, key, expr) {
  suffix <- getOption("optop.revision_suffix", "")
  if (!nzchar(suffix)) return(force(expr))
  identity <- revision_identity(cfg)
  root <- p_data("Checkpoints", suffix, stage, cfg_hash(cfg, 16L))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(root, paste0(gsub("[^A-Za-z0-9_.-]", "_", key), ".qs2"))
  if (file.exists(path)) {
    saved <- qs2::qs_read(path)
    if (!identical(saved$identity, identity)) {
      what <- c(config = !identical(saved$identity$config, identity$config),
                package = !identical(saved$identity$provenance, identity$provenance),
                code = !identical(saved$identity$code, identity$code))
      stop("Checkpoint identity mismatch (", paste(names(what)[what], collapse = ", "),
           " changed): ", path, "\n  Use a new out_suffix or remove ", root,
           call. = FALSE)
    }
    return(saved$value)
  }
  value <- force(expr)
  tmp <- tempfile("checkpoint-", tmpdir = root)
  on.exit(unlink(tmp), add = TRUE)
  qs2::qs_save(list(identity = identity, value = value), tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish checkpoint: ", path)
  value
}

# --- support resolution ----------------------------------------------------------------

#' How much of each document the harmonised support resolves, for one candidate.
#'   retained_cells         words scored individually (outside the residual bin);
#'   observed_pooled_count  scored tokens that fall in the residual bin;
#'   observed_pooled_share  the same as a share of the document's scored tokens;
#'   predicted_pooled_share fitted probability mass the candidate puts in the bin.
#' The partition is common to the candidates, so the first three do not depend
#' on K. CORPUS-level shares must pool counts (sum observed / sum L, and
#' sum L * predicted / sum L): a mean of document shares is a different quantity
#' -- see summarise_resolution().
support_resolution <- function(theta, phi, dtm, part, K) {
  part <- .densify_partition(part, dtm)
  L <- as.numeric(Matrix::rowSums(dtm))
  rare <- part$rare_mask
  stopifnot(identical(dim(rare), dim(dtm)), nrow(theta) == nrow(dtm),
            ncol(phi) == ncol(dtm))
  Tm <- as(as(dtm, "generalMatrix"), "TsparseMatrix")
  observed <- numeric(nrow(dtm))
  hit <- rare[cbind(Tm@i + 1L, Tm@j + 1L)]
  if (any(hit)) {
    a <- rowsum(Tm@x[hit], Tm@i[hit] + 1L)
    observed[as.integer(rownames(a))] <- a[, 1L]
  }
  predicted <- numeric(nrow(dtm))
  for (s in seq.int(1L, nrow(dtm), by = 256L)) {
    ii <- s:min(s + 255L, nrow(dtm))
    predicted[ii] <- rowSums((theta[ii, , drop = FALSE] %*% phi) *
                               rare[ii, , drop = FALSE])
  }
  data.table(K = K, doc_id = rownames(dtm), L = L,
             retained_cells = ncol(dtm) - rowSums(rare),
             observed_pooled_count = observed,
             observed_pooled_share = ifelse(L > 0, observed / L, NA_real_),
             predicted_pooled_share = predicted)
}

#' Corpus-level summary of support_resolution() rows, by the columns in `by`.
#' Reports BOTH the pooled token shares (the corpus quantity) and the mean of the
#' document shares (the typical document), which must not be confused.
summarise_resolution <- function(res_dt, by = "K") {
  stopifnot(is.data.table(res_dt))
  res_dt[, .(
    n_docs = .N,
    retained_cells_mean = mean(retained_cells),
    retained_cells_median = as.numeric(stats::median(retained_cells)),
    retained_cells_min = min(retained_cells),
    observed_pooled_token_share = sum(observed_pooled_count) / sum(L),
    predicted_pooled_token_share = sum(L * predicted_pooled_share) / sum(L),
    observed_pooled_doc_mean = mean(observed_pooled_share, na.rm = TRUE),
    predicted_pooled_doc_mean = mean(predicted_pooled_share)
  ), by = by]
}
