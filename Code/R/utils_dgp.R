# =============================================================================
# utils_dgp.R
# LDA data-generating process + misspecification injectors for Section 5.
#
# The vanilla DGP matches the paper (Section 5.1): phi_k ~ Dir(beta 1_W),
# theta_j ~ Dir(alpha 1_K*), N_j ~ Mult(L_j, theta_j' Phi). Counts are built in
# document chunks and stored sparse so evaluation sets of 10^4-10^5 documents
# (E2's conditional Monte Carlo) stay within memory.
#
# Injectors return generating probability rows or Theta matrices so that all
# alternatives flow through the same count sampler:
#   * contaminate_probs()     - shared "function word" distribution, doc weights
#   * burstiness (tau)        - Dirichlet-multinomial overdispersion in sampler
#   * drift_phi()             - evaluation-side topic perturbation
#   * theta_logistic_normal() - correlated topics (CTM-style)
#   * theta_mixed_alpha()     - concentration-heterogeneous document groups (E3-C)
# =============================================================================

library(Matrix)
library(data.table)

`%||dgp%` <- function(a, b) if (is.null(a)) b else a

# --- Dirichlet draws ----------------------------------------------------------

#' n draws from Dir(shape * 1_d), rows normalized; guards underflow rows.
rdirichlet_mat <- function(n, d, shape) {
  x <- matrix(rgamma(n * d, shape = shape, rate = 1), nrow = n, ncol = d)
  rs <- rowSums(x)
  bad <- which(rs < 1e-300 | !is.finite(rs))
  if (length(bad)) {
    # numerically degenerate row: all mass on one random coordinate
    x[bad, ] <- 0
    x[cbind(bad, sample.int(d, length(bad), replace = TRUE))] <- 1
    rs <- rowSums(x)
  }
  x / rs
}

# --- Document lengths ---------------------------------------------------------

#' spec: list(type = "fixed"|"poisson"|"mixture", ...); lengths floored at 10.
sample_doc_lengths <- function(spec, J, seed) {
  set.seed(seed)
  L <- switch(
    spec$type,
    fixed   = rep(as.integer(spec$L), J),
    poisson = rpois(J, spec$lambda),
    mixture = {
      n_long  <- round(J * spec$share_long)
      n_short <- J - n_long
      l <- c(rpois(n_short, spec$lambda_short), rpois(n_long, spec$lambda_long))
      grp <- c(rep("short", n_short), rep("long", n_long))
      o <- sample.int(J)                       # interleave short and long docs
      attr_l <- l[o]; attr(attr_l, "length_group") <- grp[o]
      attr_l
    },
    stop("unknown length spec type: ", spec$type)
  )
  grp <- attr(L, "length_group")
  L <- pmax(as.integer(L), 10L)
  if (!is.null(grp)) attr(L, "length_group") <- grp
  L
}

# --- Sparse multinomial sampler (chunked) --------------------------------------

#' Multinomial counts row-by-row from a dense chunk of probability rows.
#' Returns triplet components relative to chunk-local row indices.
.counts_chunk <- function(P_chunk, lengths) {
  n <- nrow(P_chunk)
  is <- vector("list", n); js <- vector("list", n); xs <- vector("list", n)
  for (r in seq_len(n)) {
    cnt <- as.integer(rmultinom(1L, size = lengths[r], prob = P_chunk[r, ]))
    nz  <- which(cnt > 0L)
    is[[r]] <- rep.int(r, length(nz)); js[[r]] <- nz; xs[[r]] <- cnt[nz]
  }
  list(i = unlist(is), j = unlist(js), x = unlist(xs))
}

#' Dirichlet-multinomial perturbation of one chunk of probability rows:
#' p~_j ~ Dir(tau * p_j). tau = Inf leaves rows unchanged (plain multinomial).
.burst_chunk <- function(P_chunk, tau) {
  if (!is.finite(tau)) return(P_chunk)
  out <- matrix(rgamma(length(P_chunk), shape = tau * P_chunk, rate = 1),
                nrow = nrow(P_chunk))
  rs <- rowSums(out)
  zero <- rs < 1e-300
  if (any(zero)) { out[zero, ] <- P_chunk[zero, , drop = FALSE]; rs[zero] <- 1 }
  out / rs
}

# --- Main corpus simulator ------------------------------------------------------

#' Simulate an LDA corpus (optionally misspecified).
#'
#' @param theta_sampler NULL for Dir(alpha); otherwise function(J, K, seed) ->
#'   matrix or list(Theta, group) for group-structured alternatives.
#' @param Phi optional fixed topic matrix (reuse the training Phi on the
#'   evaluation side; supply drift_phi(Phi, delta) output for the drift design).
#' @param contamination NULL or list(pi_stop, w_mean, w_conc).
#' @param burst_tau Dirichlet-multinomial concentration; Inf = correctly specified.
#' @param group_vocab NULL or list(doc_group [len J], block_of [len W], weight):
#'   each document mixes weight `weight` of its group's own vocabulary block into
#'   the topic mixture -> residuals concentrate by group x block (E4 alternative).
#' @return list(dtm dgCMatrix [J x W], Phi, Theta, doc_lengths, doc_group, params)
sim_lda_corpus <- function(J, W, K_true, alpha, beta, length_spec, seed,
                           Phi = NULL, theta_sampler = NULL,
                           contamination = NULL, burst_tau = Inf,
                           group_vocab = NULL,
                           chunk = 1000L, doc_prefix = "doc") {
  doc_lengths <- sample_doc_lengths(length_spec, J, seed = seed + 1L)
  doc_group   <- attr(doc_lengths, "length_group")

  set.seed(seed + 2L)
  if (is.null(Phi)) Phi <- rdirichlet_mat(K_true, W, beta)

  set.seed(seed + 3L)
  if (is.null(theta_sampler)) {
    Theta <- rdirichlet_mat(J, K_true, alpha)
  } else {
    ts <- theta_sampler(J, K_true, seed + 3L)
    if (is.list(ts)) { Theta <- ts$Theta; doc_group <- ts$group } else Theta <- ts
  }

  w_contam <- NULL; C_doc <- NULL
  if (!is.null(contamination)) {
    set.seed(seed + 4L)
    a <- contamination$w_conc * contamination$w_mean
    b <- contamination$w_conc * (1 - contamination$w_mean)
    w_contam <- if (contamination$w_mean <= 0) rep(0, J) else rbeta(J, a, b)
    if (identical(contamination$mode %||dgp% "shared", "docvary")) {
      # Doc-specific stopword mixtures: each document mixes its OWN random
      # subset of the pool with its own weights. NOT expressible as one extra
      # LDA topic (a shared pi_stop is), so this is a genuine misspecification.
      pool <- which(contamination$pi_stop > 0)
      m <- min(contamination$m %||dgp% 10L, length(pool))
      C_doc <- matrix(0, nrow = J, ncol = W)
      for (j in seq_len(J)) {
        idx <- sample(pool, m)
        wts <- rgamma(m, 1); wts <- wts / sum(wts)
        C_doc[j, idx] <- wts
      }
    }
  }

  set.seed(seed + 5L)
  starts <- seq(1L, J, by = chunk)
  trips  <- vector("list", length(starts))
  for (ci in seq_along(starts)) {
    rows <- starts[ci]:min(starts[ci] + chunk - 1L, J)
    P <- Theta[rows, , drop = FALSE] %*% Phi
    if (!is.null(w_contam)) {
      Cpart <- if (is.null(C_doc)) {
        tcrossprod(w_contam[rows], contamination$pi_stop)
      } else {
        C_doc[rows, , drop = FALSE] * w_contam[rows]
      }
      P <- P * (1 - w_contam[rows]) + Cpart
    }
    if (!is.null(group_vocab)) {
      w <- group_vocab$weight
      for (rr in seq_along(rows)) {
        g <- group_vocab$doc_group[rows[rr]]
        blk <- which(group_vocab$block_of == g)
        add <- numeric(W); add[blk] <- 1 / length(blk)
        P[rr, ] <- (1 - w) * P[rr, ] + w * add
      }
    }
    P <- .burst_chunk(P, burst_tau)
    tr <- .counts_chunk(P, doc_lengths[rows])
    tr$i <- tr$i + rows[1L] - 1L
    trips[[ci]] <- tr
  }

  dtm <- sparseMatrix(
    i = unlist(lapply(trips, `[[`, "i")),
    j = unlist(lapply(trips, `[[`, "j")),
    x = unlist(lapply(trips, `[[`, "x")),
    dims = c(J, W),
    dimnames = list(paste0(doc_prefix, "_", seq_len(J)),
                    paste0("word_", seq_len(W)))
  )

  list(
    dtm = dtm, Phi = Phi, Theta = Theta,
    doc_lengths = as.integer(doc_lengths), doc_group = doc_group,
    params = list(J = J, W = W, K_true = K_true, alpha = alpha, beta = beta,
                  length_spec = length_spec, seed = seed,
                  burst_tau = burst_tau,
                  contamination_w = if (is.null(contamination)) 0 else contamination$w_mean)
  )
}

# --- Misspecification building blocks -------------------------------------------

#' Shared "function word" distribution over n_stop random words (power-law mass).
make_stopword_dist <- function(W, n_stop = 50L, seed = 1L) {
  set.seed(seed)
  idx <- sample.int(W, n_stop)
  p <- numeric(W)
  p[idx] <- 1 / seq_len(n_stop)          # Zipf-like weights over the block
  p / sum(p)
}

#' Document-group x vocabulary-block assignment for the E4 group_vocab
#' alternative: G groups of documents, G contiguous vocabulary blocks.
make_group_vocab <- function(J, W, G = 4L, weight = 0.1, seed = 1L) {
  set.seed(seed)
  doc_group <- sample(rep_len(seq_len(G), J))
  block_of <- rep_len(seq_len(G), W)          # block g = columns with block_of==g
  list(doc_group = doc_group, block_of = block_of, weight = weight, G = G)
}

#' Evaluation-side topic drift: rows of (1-delta)*Phi + delta*Phi_fresh.
drift_phi <- function(Phi, delta, beta, seed) {
  if (delta <= 0) return(Phi)
  set.seed(seed)
  fresh <- rdirichlet_mat(nrow(Phi), ncol(Phi), beta)
  out <- (1 - delta) * Phi + delta * fresh
  out / rowSums(out)
}

#' Correlated topics: logistic-normal theta with exchangeable correlation rho.
theta_logistic_normal <- function(rho, sigma = 1) {
  function(J, K, seed) {
    set.seed(seed)
    Sig <- sigma^2 * ((1 - rho) * diag(K) + rho * matrix(1, K, K))
    Z <- matrix(rnorm(J * K), J, K) %*% chol(Sig)
    E <- exp(Z - apply(Z, 1L, max))
    E / rowSums(E)
  }
}

#' Document groups with different Dirichlet concentrations (E3 scenario C):
#' low-alpha docs are near single-topic (far from the corpus baseline),
#' high-alpha docs are mixed (close to it) -- an atypicality channel by design.
theta_mixed_alpha <- function(alphas = c(0.2, 2), share_first = 0.5) {
  function(J, K, seed) {
    set.seed(seed)
    n1 <- round(J * share_first)
    grp <- sample(c(rep("concentrated", n1), rep("diffuse", J - n1)))
    Theta <- matrix(NA_real_, J, K)
    Theta[grp == "concentrated", ] <-
      rdirichlet_mat(sum(grp == "concentrated"), K, alphas[1])
    Theta[grp == "diffuse", ] <-
      rdirichlet_mat(sum(grp == "diffuse"), K, alphas[2])
    list(Theta = Theta, group = grp)
  }
}
