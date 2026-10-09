# =============================================================================
# run_revision_diagnostics.R -- additional evaluation-only analyses of the
# revision. Cached fits only (optop.no_fit is set: a cache miss is an error).
#
# mode=e5        Support resolution of every E5 design (nested grids at the
#                configured c; the c-sensitivity arm at delta = 1): individually
#                retained cells, observed and predicted mass pooled into the
#                residual bin, and the floor-excluded document share -- reported
#                SEPARATELY, because they are different mechanisms.
#
# mode=unbinned  Protocol-matched UNBINNED comparator. For each source (MD&A;
#                seed 1 of each E1 generating configuration), candidate K,
#                protocol (reconstruction / completion) and smoothing floor f:
#                  phi_f = renormalise(max(phi, f)),  pi_f = renormalise(max(pi, f)),
#                theta folded in by the same fixed-phi ML EM on the protocol's
#                fold-in tokens under phi_f, and the SCORED tokens evaluated on
#                the full vocabulary -- no rare-word bin:
#                  D_j(model) = 2 sum_w N log(N / (L p_f)),  log score = sum_w N log p_f.
#                The unbinned deviance is finite only if every observed cell has
#                strictly positive probability; the unsmoothed WarpLDA phi is
#                mostly exact zeros, so the floor IS the convention, and it is
#                reported at 1e-12 with 1e-14 and 1e-10 as sensitivity.
#
# Usage:
#   Rscript Code/run_revision_diagnostics.R [workers] mode=e5 input=<E5 results.qs2>
#   Rscript Code/run_revision_diagnostics.R [workers] mode=unbinned
#           [sources=<E1 results.qs2>[,...]] [mdna=1] [mdna_sample_n=0]
#           [mdna_grid=10:200:10] [out_suffix=_rev2]
# =============================================================================

suppressMessages({library(data.table); library(Matrix); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
ov <- cli$overrides
SFX <- ov$out_suffix %||% "_rev2"
MODE <- ov$mode %||% "e5"
stopifnot(MODE %in% c("e5", "unbinned"))
options(optop.no_fit = !grepl("_smoke", SFX) || nzchar(Sys.getenv("OPTOP_NO_FIT")),
        optop.revision_suffix = SFX)
setup_parallel(cli$workers)
log_msg("=== revision diagnostics [%s] suffix=%s ===", MODE, SFX)

#' Seed-1 training corpus + cached fits of a simulation result object.
load_sim <- function(path, grid = NULL, extra = NULL) {
  if (!file.exists(path)) stop("input not found: ", path, call. = FALSE)
  cfg <- qs_read(path)$config
  seeds <- make_seeds(cfg$seed_base, 1L)
  tr <- sim_lda_corpus(cfg$J_train, cfg$W, cfg$K_true, cfg$alpha_DGP,
                       cfg$beta_DGP, cfg$length_spec, seed = seeds$dgp_seed,
                       doc_prefix = "tr")
  ks <- grid %||% cfg$K_grid
  fits <- get_fits_cached(tr$dtm, ks, cfg$fit_method, cfg$n_starts,
                          seeds$fit_seed_base,
                          dgp_signature(cfg, seeds, extra = extra))$models
  list(cfg = cfg, seeds = seeds, tr = tr, fits = fits)
}

# =================================== mode: e5 =====================================
if (MODE == "e5") {
  path <- ov$input
  if (is.null(path)) stop("input=<explicit E5 results.qs2> is required", call. = FALSE)
  cfg0 <- qs_read(path)$config
  z <- load_sim(path, grid = sort(unique(unlist(cfg0$grids))), extra = "e5")
  cfg <- z$cfg
  # the designs of run_E5.R: nested grids at the configured c with the package's
  # default floor (delta = c), and the c arm on the main grid at delta = 1
  designs <- c(
    lapply(names(cfg$grids), function(nm)
      list(name = paste0("grid_", nm), ks = cfg$grids[[nm]], c = cfg$c,
           min_null = NULL)),
    lapply(cfg$c_grid, function(cc)
      list(name = paste0("c_", num2tag(cc)), ks = cfg$grids[[1L]], c = cc,
           min_null = 1)))
  rows <- lapply(designs, function(d)
    revision_checkpoint("E5_resolution", cfg, d$name, {
      fl <- z$fits[as.character(d$ks)]
      sc <- score_insample(fl, z$tr$dtm, d$c, c("dev", "chisq"),
                           min_null = d$min_null)
      rr <- rbindlist(lapply(names(fl), function(k)
        support_resolution(theta_from_fit(fl[[k]]), phi_from_fit(fl[[k]]),
                           z$tr$dtm, sc$partition, as.integer(k))))
      log_msg("E5 design '%s' (c = %s, |K| = %d) resolved", d$name, d$c, length(fl))
      list(resolution = rr[, `:=`(design = d$name, c = d$c)],
           floor = copy(sc$summary)[, `:=`(design = d$name, c = d$c)],
           minbin = data.table(design = d$name, c = d$c,
                               share = sc$minbin_report$share,
                               excluded_mass = sc$minbin_report$excluded_mass))
    }))
  out <- list(resolution = rbindlist(lapply(rows, `[[`, "resolution")),
              floor = rbindlist(lapply(rows, `[[`, "floor")),
              minbin = rbindlist(lapply(rows, `[[`, "minbin")),
              config = cfg, source = path)
  out$summary <- summarise_resolution(out$resolution, by = c("design", "c", "K"))
  f <- p_data("E5", paste0("e5_resolution", SFX, ".qs2"))
  cache_put(out, f, cfg)
  fwrite(out$summary, p_results("csv", paste0("e5_resolution_summary", SFX, ".csv")))
  print(out$summary)
  log_msg("saved %s", basename(f))
}

# ================================= mode: unbinned ==================================
if (MODE == "unbinned") {
  FLOORS <- c(1e-14, 1e-12, 1e-10)
  floor_renorm <- function(p, f) { p <- pmax(p, f); p / sum(p) }

  #' Fixed-phi ML EM in document blocks (documents are independent given phi, so
  #' blocking is exact up to the block-wise stopping rule) -- bounds the dense
  #' J x W work matrix on the MD&A corpus.
  foldin_blocked <- function(phi, dtm, block = 1000L) {
    do.call(rbind, lapply(split(seq_len(nrow(dtm)),
                                ceiling(seq_len(nrow(dtm)) / block)),
                          function(ii) foldin_theta_em(phi, dtm[ii, , drop = FALSE])))
  }

  score_unbinned <- function(src, nm, k, protocol, ff) {
    fold <- if (protocol == "rec") src$ev else src$spl$foldin
    score <- if (protocol == "rec") src$ev else src$spl$score
    phi0 <- phi_from_fit(src$fits[[k]])
    ph <- pmax(phi0, ff); ph <- ph / rowSums(ph)
    pf <- floor_renorm(src$pi, ff)
    th <- foldin_blocked(ph, fold)
    L <- as.numeric(Matrix::rowSums(score))
    dm <- dn <- ll <- numeric(nrow(score)); n_zero <- 0; tok_zero <- 0
    for (s in seq.int(1L, nrow(score), by = 256L)) {
      ii <- s:min(s + 255L, nrow(score))
      Tm <- as(as(score[ii, , drop = FALSE], "generalMatrix"), "TsparseMatrix")
      if (!length(Tm@x)) next
      cell <- cbind(Tm@i + 1L, Tm@j + 1L)
      logp <- log((th[ii, , drop = FALSE] %*% ph)[cell])
      # observed cells on which the floor BINDS: the unsmoothed predictor puts
      # less than f there (exact zeros included), so their score is set by f
      zero <- (th[ii, , drop = FALSE] %*% phi0)[cell] < ff
      n_zero <- n_zero + sum(zero); tok_zero <- tok_zero + sum(Tm@x[zero])
      Li <- L[ii][Tm@i + 1L]
      a <- rowsum(cbind(dm = 2 * Tm@x * (log(Tm@x / Li) - logp),
                        dn = 2 * Tm@x * (log(Tm@x / Li) - log(pf[Tm@j + 1L])),
                        ll = Tm@x * logp), Tm@i + 1L)
      at <- ii[as.integer(rownames(a))]
      dm[at] <- a[, "dm"]; dn[at] <- a[, "dn"]; ll[at] <- a[, "ll"]
    }
    scoring_log("unbinned %s K=%s %s floor=%g done", nm, k, protocol, ff)
    list(doc = data.table(source = nm, K = as.integer(k), protocol = protocol,
                          floor = ff, doc_id = rownames(score), L = L,
                          d_model = dm, d_null = dn, log_score = ll),
         zero = data.table(source = nm, K = as.integer(k), protocol = protocol,
                           floor = ff, n_floor_cells = n_zero, tokens_floor = tok_zero,
                           tokens = sum(L)))
  }

  sources <- list()
  if ((ov$mdna %||% "1") == "1") {
    n_s <- as.integer(ov$mdna_sample_n %||% "0")
    prep <- qs_read(p_data("MDNA", sprintf("mdna_prep_2015_2016%s.qs2",
                                           if (n_s > 0L) sprintf("_n%d", n_s) else "")))
    tr <- prep$dtm_train
    sig <- list(corpus = "mdna_item7_pooled", y1 = 2015L, y2 = 2016L,
                dtm_hash = cfg_hash(list(dim(tr), Matrix::rowSums(tr)[1:20],
                                         colnames(tr)[1:50])))
    grid <- .parse_override_value(ov$mdna_grid %||% "10:200:10")
    sources$MDNA <- list(
      ev = prep$dtm_ev, pi = as.numeric(OpTop::optop_make_baseline(tr)$pi_glob),
      fits = get_fits_cached(tr, grid, "WarpLDA", 3L, 1970L, sig)$models,
      spl = split_tokens_binomial(prep$dtm_ev, 0.5, seed = 1970L))
    rm(prep, tr)
  }
  src_paths <- if (!is.null(ov$sources)) strsplit(ov$sources, ",", fixed = TRUE)[[1L]] else
    p_data("E1", sprintf("e1_results_E1_full_%s.qs2", c(
      "Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01",
      "Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01")))
  for (path in src_paths) {
    z <- load_sim(path)
    ev <- sim_lda_corpus(z$cfg$J_eval, z$cfg$W, z$cfg$K_true, z$cfg$alpha_DGP,
                         z$cfg$beta_DGP, z$cfg$length_spec,
                         seed = z$seeds$dgp_seed + 7L, Phi = z$tr$Phi,
                         doc_prefix = "ev")$dtm          # as run_E1.R, replicate 1
    sources[[sprintf("E1_Kstar%d_W%d", z$cfg$K_true, z$cfg$W)]] <- list(
      ev = ev, pi = as.numeric(OpTop::optop_make_baseline(z$tr$dtm)$pi_glob),
      fits = z$fits,
      spl = split_tokens_binomial(ev, z$cfg$completion_prop,
                                  seed = z$seeds$split_seed))
  }

  res <- list()
  for (nm in names(sources)) {
    src <- sources[[nm]]
    tasks <- CJ(ff = FLOORS, k = names(src$fits), protocol = c("rec", "com"),
                sorted = FALSE)
    ident <- list(source = nm, eval = digest::digest(src$ev, algo = "xxhash64"),
                  grid = names(src$fits))
    log_msg("unbinned: %s -- %d tasks", nm, nrow(tasks))
    res[[nm]] <- future_lapply(seq_len(nrow(tasks)), function(i)
      revision_checkpoint("unbinned", ident,
                          sprintf("K%s_%s_f%g", tasks$k[i], tasks$protocol[i], tasks$ff[i]),
                          score_unbinned(src, nm, tasks$k[i], tasks$protocol[i],
                                         tasks$ff[i])),
      future.seed = NULL)
  }
  flat <- unlist(res, recursive = FALSE)
  docs <- rbindlist(lapply(flat, `[[`, "doc"))
  zero <- rbindlist(lapply(flat, `[[`, "zero"))
  # delta = 1 on the unbinned null discrepancy, as everywhere else in the paper
  docs[, r2_doc := ifelse(d_null >= 1, 1 - d_model / d_null, NA_real_)]
  summ <- docs[, .(r2_micro = 1 - sum(d_model[!is.na(r2_doc)]) / sum(d_null[!is.na(r2_doc)]),
                   r2_macro = mean(r2_doc, na.rm = TRUE),
                   log_score_per_token = sum(log_score) / sum(L),
                   perplexity = exp(-sum(log_score) / sum(L)),
                   excluded_share = mean(is.na(r2_doc)), n_docs = .N),
               by = c("source", "protocol", "floor", "K")]
  summ <- merge(summ, zero[, .(source, protocol, floor, K, n_floor_cells,
                               token_share_floor = tokens_floor / tokens)],
                by = c("source", "protocol", "floor", "K"))
  out <- list(doc = docs, summary = summ, floors = FLOORS,
              convention = paste("phi and the baseline floored at f and renormalised;",
                                 "theta refolded under the floored phi; delta = 1"),
              scope = "MD&A full grid; seed 1 of each E1 generating configuration")
  f <- p_data("MDNA", paste0("unbinned_comparator", SFX, ".qs2"))
  cache_put(out, f, list(experiment = "unbinned", floors = FLOORS))
  fwrite(summ, p_results("csv", paste0("unbinned_comparator_summary", SFX, ".csv")))
  print(summ[floor == 1e-12, .(K_best_log_score = K[which.max(log_score_per_token)],
                               K_best_micro = K[which.max(r2_micro)]),
             by = c("source", "protocol")])
  log_msg("saved %s", basename(f))
}
log_msg("=== revision diagnostics [%s] complete ===", MODE)
