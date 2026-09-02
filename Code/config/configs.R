# =============================================================================
# configs.R
# One config constructor per experiment, three profiles each:
#   smoke - minutes-scale end-to-end exercise of every module (tiny corpus);
#   pilot - <= 1 h wall on ~11 workers; every figure/table gets real data;
#   full  - the paper-scale design (two overnight runs).
# All parameters live here; run tags/caches derive from these lists.
# =============================================================================

.base_dgp <- function() {
  list(
    K_true = 10L, J_train = 1000L, J_eval = 500L, W = 10000L,
    alpha_DGP = 0.5, beta_DGP = 0.01,     # data-generating Dirichlet priors
    alpha = NA_real_, beta = NA_real_,    # model/fit priors (NA = engine default)
    length_spec = list(type = "fixed", L = 1000L),
    K_grid = 2:20, c = 1, n_starts = 3L,
    fit_method = "VEM",              # VEM | Gibbs | WarpLDA (text2vec)
    metrics = c("dev", "chisq", "se"),
    completion_prop = 0.5, sel_alpha = 0.05, eps_grid = c(0.01, 0.005),
    seed_base = 1970000L
  )
}

.smoke_dgp <- list(J_train = 120L, J_eval = 60L, W = 800L, K_true = 4L,
                   K_grid = 2:6, n_starts = 1L)

cfg_E1 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E1", profile = profile), .base_dgp())
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(S = 2L)),
    pilot = list(S = 3L, n_starts = 2L, K_grid = 2:15),
    full  = list(S = 50L)
  )
  cfg <- utils::modifyList(cfg, extra)
  # DEPRECATED, no effect: the Gibbs robustness arm was removed in the 4-study
  # consolidation. Kept so existing `gibbs_S=0` CLI overrides don't error.
  cfg$gibbs_S <- 0L
  cfg
}

cfg_E2 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E2", profile = profile), .base_dgp())
  cfg$metrics <- "dev"
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(S_train = 1L, R_eval = 10L, J_ev_grid = c(40L, 80L),
                               J_truth = 400L, K_test = c(3L, 4L, 5L),
                               n_starts = 1L)),
    pilot = list(S_train = 2L, R_eval = 60L, J_ev_grid = c(100L, 400L),
                 J_truth = 4000L, K_test = c(5L, 10L, 15L), n_starts = 2L,
                 K_grid = 2:15),
    # R = 300 gives MC-SE ~1.3% on coverage; J_ev = 1000 dropped (CLT already
    # demonstrated at 500); calibrated to ~50 core-h on the reference machine.
    full  = list(S_train = 10L, R_eval = 300L,
                 J_ev_grid = c(100L, 250L, 500L),
                 J_truth = 20000L, K_test = c(5L, 10L, 15L))
  )
  utils::modifyList(cfg, extra)
}

cfg_E3 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E3", profile = profile), .base_dgp())
  cfg$scenarios <- list(
    A = list(label = "A: Heterogeneous lengths",
             length_spec = list(type = "mixture", share_long = 0.2,
                                lambda_short = 500, lambda_long = 5000),
             theta = "dirichlet"),
    B = list(label = "B: Homogeneous lengths",
             length_spec = list(type = "poisson", lambda = 1000),
             theta = "dirichlet"),
    C = list(label = "C: Atypicality (mixed concentration)",
             length_spec = list(type = "fixed", L = 1000L),
             theta = "mixed_alpha", alphas = c(0.2, 2), share_first = 0.5)
  )
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(S = 1L, n_starts = 1L)),
    pilot = list(S = 2L, n_starts = 1L, K_grid = 2:15),
    # single-start: E3 estimates gaps/decompositions, not K selection, and
    # scenario-A corpora (5000-token docs) fit ~3x slower; multi-start
    # robustness is documented in E1. ~85 core-h at S = 15.
    full  = list(S = 15L, n_starts = 1L)
  )
  utils::modifyList(cfg, extra)
}

cfg_E4 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E4", profile = profile), .base_dgp())
  cfg$metrics <- "dev"
  cfg$B_strata <- 5L; cfg$S_strata <- 5L; cfg$min_docfreq <- 5L
  cfg$n_stop <- 50L
  cfg$n_groups <- 4L                 # document groups for the group_vocab alt
  full_alts <- list(
    contamination      = list(strengths = c(0.05, 0.10, 0.20)),  # docvary, both sides
    contamination_eval = list(strengths = c(0.05, 0.10, 0.20)),  # shared, eval-only
    burstiness         = list(strengths = c(2000, 500, 100)),
    drift              = list(strengths = c(0.05, 0.10, 0.20)),
    ctm                = list(strengths = c(0.3, 0.6, 0.9)),
    group_vocab        = list(strengths = c(0.05, 0.10, 0.20))   # doc-group vocab (§5.2)
  )
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(
      S_train = 1L, R_null = 10L, R_power = 5L, J_ev = 60L, J_center = 300L,
      K_test = c(3L, 4L), n_starts = 1L,
      alternatives = list(group_vocab = list(strengths = 0.2)))),
    pilot = list(S_train = 2L, R_null = 60L, R_power = 40L, J_ev = 500L,
                 J_center = 2000L,
                 K_test = c(5L, 10L, 15L), n_starts = 2L, K_grid = 2:15,
                 alternatives = list(
                   contamination = list(strengths = c(0.05, 0.20)),
                   contamination_eval = list(strengths = 0.10),
                   group_vocab = list(strengths = c(0.05, 0.20)))),
    full  = list(S_train = 10L, R_null = 500L, R_power = 200L, J_ev = 500L,
                 J_center = 10000L,
                 K_test = c(5L, 10L, 15L), alternatives = full_alts)
  )
  utils::modifyList(cfg, extra)
}

cfg_E5 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E5", profile = profile), .base_dgp())
  cfg$metrics <- c("dev", "chisq")
  cfg$c_grid <- c(1, 5)
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(grids = list(g1 = 2:5, g2 = 2:6),
                               K_grid = 2:6, n_starts = 1L)),
    pilot = list(grids = list(g2_15 = 2:15, g2_22 = 2:22),
                 K_grid = 2:22, n_starts = 1L),
    # nested grids CONTAINING the paper design point K* = 40 (WarpLDA makes
    # the large-K fits cheap; the old 2:20/2:50 grids sat entirely below K*)
    full  = list(grids = list(g10_100 = seq(10L, 100L, 10L),
                              g10_160 = seq(10L, 160L, 10L)),
                 K_grid = seq(10L, 160L, 10L), n_starts = 1L)
  )
  utils::modifyList(cfg, extra)
}

cfg_E6 <- function(profile = c("pilot", "smoke", "full")) {
  profile <- match.arg(profile)
  cfg <- c(list(experiment = "E6", profile = profile), .base_dgp())
  cfg$metrics <- "dev"                     # word-level deviance is the primary index
  cfg$word_filter <- list(min_docfreq = 5L, min_expected = 5)
  # Two diagnostic regimes for the w-Micro/w-Macro divergence (§3.8.3):
  #   W1 correctly specified -> model fits common words well but the (rare,
  #      estimation-hard) vocabulary poorly, so w-Micro >> w-Macro;
  #   W2 high-frequency stopword contamination -> those non-topic common words
  #      are mispredicted (r2_w ~ 0) and, carrying large baseline-discrepancy
  #      weight, pull w-Micro DOWN toward w-Macro; they are the worst-fit
  #      common words the diagnostic recovers.
  cfg$scenarios <- list(
    W1 = list(label = "W1: Correctly specified", contam = FALSE),
    W2 = list(label = "W2: High-frequency stopword contamination",
              contam = TRUE, w_mean = 0.05)
  )
  extra <- switch(profile,
    smoke = c(.smoke_dgp, list(S = 1L, n_starts = 1L, n_stop = 40L)),
    pilot = list(S = 2L, n_starts = 1L, K_grid = 2:15, n_stop = 60L),
    full  = list(S = 15L, n_stop = 80L)
  )
  utils::modifyList(cfg, extra)
}

get_config <- function(experiment, profile) {
  switch(experiment,
    E1 = cfg_E1(profile), E2 = cfg_E2(profile), E3 = cfg_E3(profile),
    E4 = cfg_E4(profile), E5 = cfg_E5(profile), E6 = cfg_E6(profile),
    # E1 robustness variants (todo.txt): higher K* / sparser doc-topic prior.
    # MANUAL ONLY -- not in make_all.R's full orchestration since the 4-study
    # consolidation. Replications and starts capped because the 5:40 grid is
    # expensive (K = 40 fits ~19 min each; ~90 core-h at S = 10).
    E1b = utils::modifyList(cfg_E1(profile), list(
      experiment = "E1b", K_true = 20L, K_grid = 5:40,
      S = min(cfg_E1(profile)$S, 10L), n_starts = 2L)),
    E1c = utils::modifyList(cfg_E1(profile), list(
      experiment = "E1c", alpha_DGP = 0.1,
      S = min(cfg_E1(profile)$S, 15L), n_starts = 2L)),
    stop("unknown experiment: ", experiment)
  )
}

#' Parse a CLI override value:
#'   "2:30"               -> 2:30            (contiguous range)
#'   "5:100:5"            -> seq(5, 100, 5)  (stepped range: from:to:by)
#'   "10,20,40,50,70,100" -> that vector     (arbitrary set)
#'   plain numbers        -> numeric/integer; anything else -> string.
.parse_override_value <- function(x) {
  if (grepl("^-?[0-9]+:-?[0-9]+:[0-9]+$", x)) {
    p <- as.integer(strsplit(x, ":", fixed = TRUE)[[1L]])
    return(seq(p[1L], p[2L], by = p[3L]))
  }
  if (grepl("^-?[0-9]+:[0-9]+$", x)) {
    p <- as.integer(strsplit(x, ":", fixed = TRUE)[[1L]])
    return(p[1L]:p[2L])
  }
  if (grepl("^-?[0-9]+(,-?[0-9]+)+$", x)) {
    return(as.integer(strsplit(x, ",", fixed = TRUE)[[1L]]))
  }
  if (grepl("^-?[0-9]*\\.?[0-9]+$", x)) {
    v <- as.numeric(x)
    return(if (!grepl(".", x, fixed = TRUE)) as.integer(v) else v)
  }
  x
}

#' Apply CLI `key=value` overrides to a profile config. The shorthand `L` sets a
#' fixed document length. Note: E3's scenarios carry their own scenario-specific
#' length_specs, which `L` does not touch.
#'
#' `strict = TRUE` (default, for direct single-stage calls) errors on a key not
#' in the config. `strict = FALSE` (used by the experiment drivers and
#' make_figures/make_tables, which receive one shared override set spanning many
#' experiments) skips such keys with a one-line notice — so a shared set like
#' `S=3 gibbs_S=0` applies to the experiments that use those keys and is ignored
#' elsewhere. make_all.R's startup guard still fails fast on true typos.
apply_overrides <- function(cfg, overrides, strict = TRUE) {
  if (!length(overrides)) return(cfg)
  skipped <- character(0)
  for (key in names(overrides)) {
    val <- .parse_override_value(overrides[[key]])
    if (key == "L") {
      cfg$length_spec <- list(type = "fixed", L = as.integer(val))
      next
    }
    if (!key %in% names(cfg)) {
      if (strict) {
        stop(sprintf("unknown override '%s'.\nValid keys: L, %s",
                     key, paste(sort(names(cfg)), collapse = ", ")),
             call. = FALSE)
      }
      skipped <- c(skipped, key)
      next
    }
    cfg[[key]] <- val
  }
  if (length(skipped)) {
    who <- if (!is.null(cfg$experiment)) cfg$experiment else "this config"
    message(sprintf("apply_overrides: ignoring %d override(s) not used by %s: %s",
                    length(skipped), who, paste(skipped, collapse = ", ")))
  }
  # alpha/beta on the CLI are the FITTED-MODEL priors -- easy to confuse with
  # the data-generating priors, which are alpha_DGP/beta_DGP.
  fp <- intersect(names(overrides), c("alpha", "beta"))
  fp <- fp[vapply(fp, function(k) is.numeric(cfg[[k]]) && !is.na(cfg[[k]]),
                  logical(1L))]
  if (length(fp)) {
    message(sprintf(paste(
      "NOTE: %s set the FITTED-MODEL prior(s) (WarpLDA doc_topic_prior/",
      "topic_word_prior, Gibbs alpha/delta, VEM alpha). The data-generating",
      "priors are alpha_DGP=%s, beta_DGP=%s."),
      paste(fp, collapse = "/"), cfg$alpha_DGP, cfg$beta_DGP))
  }
  for (k in c("J_train", "J_eval", "W", "K_true", "n_starts",
              "alpha_DGP", "beta_DGP")) {
    if (!is.null(cfg[[k]]) && any(cfg[[k]] <= 0)) {
      stop(sprintf("override leaves %s <= 0", k), call. = FALSE)
    }
  }
  for (k in c("alpha", "beta")) {          # model priors: NA = engine default
    if (!is.null(cfg[[k]]) && !is.na(cfg[[k]]) && cfg[[k]] <= 0) {
      stop(sprintf("override leaves model %s <= 0", k), call. = FALSE)
    }
  }
  if (!is.null(cfg$K_true) && !is.null(cfg$K_grid) &&
      !(cfg$K_true %in% cfg$K_grid)) {
    warning(sprintf("K_true = %d is not inside K_grid [%d..%d]",
                    cfg$K_true, min(cfg$K_grid), max(cfg$K_grid)),
            call. = FALSE)
  }
  cfg
}

#' Signature identifying a simulated corpus (drives the shared fit cache).
dgp_signature <- function(cfg, seeds, scenario = "base",
                          extra = NULL) {
  sig <- list(K_true = cfg$K_true, J = cfg$J_train, W = cfg$W,
       alpha = cfg$alpha_DGP, beta = cfg$beta_DGP, length_spec = cfg$length_spec,
       dgp_seed = seeds$dgp_seed, scenario = scenario, extra = extra)
  # Model/fit priors are fit-side but key the SAME per-(corpus,K) fit cache;
  # fold them in ONLY when set so unset runs keep their existing cache hashes.
  # prefit_pool()/.fit_one_k() read these back to fit with the requested prior.
  if (!is.na(cfg$alpha)) sig$fit_alpha <- cfg$alpha
  if (!is.na(cfg$beta))  sig$fit_beta  <- cfg$beta
  sig
}
