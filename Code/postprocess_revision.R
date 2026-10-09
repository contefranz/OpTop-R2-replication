# =============================================================================
# postprocess_revision.R -- CACHE-ONLY recomputations for the Sept-2026
# revision. Reads saved result objects; touches no fit, folds in nothing,
# overwrites nothing (every output carries the suffix below).
#
#   2a  MD&A selections under three rules (pointwise adjacent [exploratory],
#       simultaneous adjacent [local], total gain [primary]) x protocol x eps
#       x family, + the full gain profile with all three upper bounds
#   2b  firm-cluster standard errors; subset of evaluation firms absent from
#       training
#   2c  discrepancy-floor (delta) sensitivity at FIXED support threshold c
#   2d  the same rules on the c = 0.5 / c = 2 / truncated-grid caches
#   2e  residual mass displaced between training-frequency groups
#   2f  simulations (baseline + second DGP): simultaneous adjacent selection
#       from cached gains, exact total-gain for the one replicate with
#       document scores, a lower bound on the total-gain selection for all
#       seeds, and the protocol-MATCHED likelihood optimum
#   2g  coverage experiment: calibration against a NOISY reference target
#
# Every step ends in a gate against values computed independently during the
# audit; a failed gate stops the script.
#
# Usage:  Rscript Code/postprocess_revision.R [out_suffix=_rev1]
# =============================================================================

suppressMessages({library(data.table); library(qs2); library(Matrix)})
source(here::here("Code", "R", "source_all.R"))
options(optop.no_fit = TRUE)

.args <- commandArgs(trailingOnly = TRUE)
SFX <- sub("^out_suffix=", "", grep("^out_suffix=", .args, value = TRUE))
if (!length(SFX)) SFX <- "_rev1"
ALPHA <- 0.05; EPS <- c(0.01, 0.005)
gate <- function(ok, fmt, ...) {
  msg <- sprintf(fmt, ...)
  if (!isTRUE(ok)) stop("GATE FAILED: ", msg, call. = FALSE)
  log_msg("GATE ok: %s", msg)
}
soft <- function(ok, fmt, ...) log_msg("%s %s", if (isTRUE(ok)) "CHECK ok:" else
  "CHECK DIFFERS (not fatal):", sprintf(fmt, ...))
near <- function(a, b, tol = 5e-9) isTRUE(all(abs(a - b) <= tol))
emit <- function(dt, name) {
  fwrite(dt, p_results("csv", sprintf("%s%s.csv", name, SFX)))
  log_msg("csv %s%s.csv (%d rows)", name, SFX, nrow(dt))
}
# share selecting exactly K*: a replicate with NO selection counts as a miss
hit <- function(k, k_true) !is.na(k) & k == k_true
kh <- function(dt, rule_v, eps_v, metric_v = "dev")
  dt[rule == rule_v & eps == eps_v & metric == metric_v, K_hat]
OUT <- list()

# ------------------------------ MD&A inputs -----------------------------------
x <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016.qs2"))
prep <- qs_read(p_data("MDNA", "mdna_prep_2015_2016.qs2"))
cik <- data.table(doc_id = prep$dv_ev$doc_id,
                  cluster_id = as.character(prep$dv_ev$cik))
docs <- list(com = x$doc_com, rec = x$doc_rec)

rules_for <- function(doc_list, cluster = NULL)
  rbindlist(lapply(names(doc_list), function(pn)
    select_k_all_rules(doc_list[[pn]], EPS, ALPHA, cluster)[, protocol := pn]))

# ================================== 2a ==========================================
sel <- rules_for(docs)
OUT$mdna_sel <- sel
for (pn in c("com", "rec")) {
  s <- sel[protocol == pn]
  gate(identical(kh(s, "adjacent_pointwise", .01), if (pn == "com") 50L else 80L) &&
       identical(kh(s, "adjacent_simultaneous", .01), 80L) &&
       identical(kh(s, "total_gain", .01), 180L) &&
       identical(kh(s, "adjacent_pointwise", .005), 120L) &&
       identical(kh(s, "adjacent_simultaneous", .005), 120L) &&
       identical(kh(s, "total_gain", .005), 190L),
       "2a MD&A Deviance %s: pointwise/sim.-adjacent/total-gain = %s/80/180 (.01), 120/120/190 (.005)",
       pn, if (pn == "com") "50" else "80")
}
# gain profile: adjacent steps with both bounds + the total-gain bound per K
prof <- rbindlist(lapply(names(docs), function(pn) {
  adj <- paired_gains_all(docs[[pn]][metric == "dev"], ALPHA, NULL, "adjacent")
  allp <- paired_gains_all(docs[[pn]][metric == "dev"], ALPHA, NULL, "all_pairs")
  tg <- allp[, .(total_gain_max = max(delta_mean), total_gain_ub = max(ub_simul),
                 K_argmax = K_to[which.max(delta_mean)], n_larger = .N), by = K]
  merge(adj[, .(K, K_to, delta_mean, se, ub_pointwise, ub_adjacent_simul = ub_simul,
                z_adjacent = z_crit)], tg, by = "K")[, protocol := pn]
}))
p5060 <- prof[protocol == "com" & K == 50]
allp_com <- paired_gains_all(x$doc_com[metric == "dev"], ALPHA, NULL, "all_pairs")
gate(near(p5060$delta_mean, 0.008109138) && near(p5060$se, 0.0009272185, 5e-11) &&
     near(p5060$ub_pointwise, 0.009634277) && near(p5060$ub_adjacent_simul, 0.01069651, 5e-9) &&
     near(allp_com[K == 50 & K_to == 200, delta_mean], 0.097829152),
     "2a completion 50->60: gain %.9f se %.10f UB %.9f Bonf.UB %.8f; 50->200 gain %.9f",
     p5060$delta_mean, p5060$se, p5060$ub_pointwise, p5060$ub_adjacent_simul,
     allp_com[K == 50 & K_to == 200, delta_mean])
OUT$mdna_profile <- prof
emit(sel[, .(protocol, metric, rule, eps, K_hat, certified, n_comparisons, max_ub,
             M, se_type)], "mdna_selection_rules")
emit(prof, "mdna_gain_profile")

# ================================== 2b ==========================================
n_firms <- uniqueN(prep$dv_ev$cik)
n_two <- sum(table(prep$dv_ev$cik) == 2L)
in_tr <- prep$dv_ev$cik %in% prep$dv_train$cik
gate(n_firms == 3059L && n_two == 407L && sum(in_tr) == 1963L,
     "2b firms: %d evaluation firms, %d with two filings, %d evaluation docs from firms also in training",
     n_firms, n_two, sum(in_tr))
sel_cl <- rules_for(docs, cluster = cik)
adj_cl <- paired_gains_all(x$doc_com[metric == "dev"], ALPHA, cik, "adjacent")
gate(near(adj_cl[K == 50, se], 0.0009863081, 5e-11),
     "2b firm-cluster SE of the completion 50->60 gain = %.10f (iid %.10f, +%.1f%%; %d clusters)",
     adj_cl[K == 50, se], p5060$se, 100 * (adj_cl[K == 50, se] / p5060$se - 1),
     adj_cl$n_clusters[1L])
same_sel <- merge(sel[, .(protocol, metric, rule, eps, K_iid = K_hat)],
                  sel_cl[, .(protocol, metric, rule, eps, K_cluster = K_hat)],
                  by = c("protocol", "metric", "rule", "eps"))
gate(same_sel[metric == "dev", all(K_iid == K_cluster)],
     "2b firm-cluster SEs leave every full-sample Deviance selection unchanged")
# evaluation firms absent from training (no new fit; changes the population)
abs_ids <- prep$dv_ev$doc_id[!in_tr]
docs_abs <- lapply(docs, function(d) d[doc_id %in% abs_ids])
sel_abs <- rules_for(docs_abs, cluster = cik)
n_abs <- vapply(docs_abs, function(d) d[metric == "dev" & K == 50 & !is.na(r2_doc), .N], 0L)
gate(length(abs_ids) == 1503L, "2b absent-firm subset: %d evaluation documents (retained: rec %d, com %d)",
     length(abs_ids), n_abs[["rec"]], n_abs[["com"]])
aud <- data.table(protocol = rep(c("rec", "com"), each = 2), eps = rep(c(.01, .005), 2),
                  pw = c(50L, 120L, 50L, 120L), sa = c(50L, 150L, 50L, 140L),
                  tg = c(190L, 190L, 190L, NA))
for (i in seq_len(nrow(aud))) {
  s <- sel_abs[protocol == aud$protocol[i]]
  got <- c(kh(s, "adjacent_pointwise", aud$eps[i]), kh(s, "adjacent_simultaneous", aud$eps[i]),
           kh(s, "total_gain", aud$eps[i]))
  soft(identical(got, c(aud$pw[i], aud$sa[i], aud$tg[i])),
       "2b absent-firm %s eps=%.3f: %s (audit %s)", aud$protocol[i], aud$eps[i],
       paste(got, collapse = "/"), paste(c(aud$pw[i], aud$sa[i], aud$tg[i]), collapse = "/"))
}
OUT$mdna_sel_cluster <- sel_cl; OUT$mdna_sel_absent <- sel_abs
emit(rbindlist(list(
  sel[metric == "dev"][, sample := "full sample, document-level SE"],
  sel_cl[metric == "dev"][, sample := "full sample, firm-cluster SE"],
  sel_abs[metric == "dev"][, sample := "firms absent from training, firm-cluster SE"]),
  use.names = TRUE)[, .(sample, protocol, rule, eps, K_hat, certified, n_comparisons,
                        max_ub, M)], "mdna_selection_firm_sensitivity")

# ================================== 2c ==========================================
# delta varied with c held at 1: the cache keeps raw d_model / d_null for every
# document, so the retained set can be rebuilt without rescoring.
refloor <- function(d, delta) {
  d <- copy(d)
  keep <- if (delta > 0) d$d_null >= delta else d$d_null > 0
  d[, r2_doc := ifelse(keep, 1 - d_model / d_null, NA_real_)]
  d
}
DELTAS <- c(0, 0.1, 0.5, 1, 2, 5)
dl <- rbindlist(lapply(DELTAS, function(dv) rbindlist(lapply(names(docs), function(pn) {
  d <- refloor(docs[[pn]][metric == "dev"], dv)
  k50 <- d[K == 50 & !is.na(r2_doc)]
  s <- select_k_all_rules(d, EPS, ALPHA)
  data.table(delta = dv, protocol = pn, retained = nrow(k50),
             excluded = d[K == 50, sum(is.na(r2_doc))],
             macro_K50 = mean(k50$r2_doc),
             micro_K50 = 1 - sum(k50$d_model) / sum(k50$d_null),
             s[, .(rule, eps, K_hat)])
}))))
dl_w <- dcast(dl, delta + protocol + retained + excluded + macro_K50 + micro_K50 ~ rule + eps,
              value.var = "K_hat")
g2c <- dl_w[protocol == "com"][order(delta)]
gate(identical(g2c$retained, c(3452L, 3449L, 3449L, 3447L, 3447L, 3441L)) &&
     near(g2c$macro_K50, c(0.273545, 0.438809, 0.438809, 0.439823, 0.439823, 0.440560), 5e-7),
     "2c delta sensitivity (completion, K=50): retained %s; Macro %s",
     paste(g2c$retained, collapse = "/"), paste(sprintf("%.4f", g2c$macro_K50), collapse = "/"))
gate(near(dl_w[delta == 1 & protocol == "com", macro_K50],
          x$summary[metric == "dev" & eval == "com" & K == 50, r2_macro], 1e-12),
     "2c delta = 1 reproduces the cached completion Macro index at K=50 (%.6f)",
     x$summary[metric == "dev" & eval == "com" & K == 50, r2_macro])
OUT$mdna_delta <- dl_w
emit(dl_w, "mdna_delta_sensitivity")

# ================================== 2d ==========================================
# These runs refined the grid around their coarse-pass selection, so their
# supports are harmonised over MORE models than the baseline's 20: the rows
# below are results on those refined supports, restricted to the coarse grid
# points for the comparison -- not a like-for-like change of c alone.
sens <- rbindlist(lapply(c(c05 = "_c05", c2 = "_c2", g100 = "_g100"), function(sf) {
  y <- qs_read(p_data("MDNA", sprintf("mdna_results_MDNA_2015_2016%s.qs2", sf)))
  coarse <- y$K_all[y$K_all %% 10L == 0L]
  dd <- list(com = y$doc_com[K %in% coarse], rec = y$doc_rec[K %in% coarse])
  s <- rules_for(dd)[metric == "dev"]
  lvl <- rbindlist(lapply(c(ins = "ins", rec = "rec", com = "com"), function(mn) {
    a <- y$summary[metric == "dev" & eval == mn & K %in% coarse, .(K, r2 = r2_macro)]
    b <- x$summary[metric == "dev" & eval == mn, .(K, r0 = r2_macro)]
    merge(a, b, by = "K")[, .(max_abs = max(abs(r2 - r0)))]
  }))
  s[, `:=`(c_part = y$params$c_part, n_models_in_support = length(y$K_all),
           n_coarse = length(coarse), max_abs_dR2 = max(lvl$max_abs))]
  s
}), idcol = "design")
for (dn in c("c05", "c2", "g100")) {
  s <- sens[design == dn]
  pw <- c(s[protocol == "com" & rule == "adjacent_pointwise" & eps == .01, K_hat],
          s[protocol == "rec" & rule == "adjacent_pointwise" & eps == .01, K_hat],
          s[protocol == "com" & rule == "adjacent_pointwise" & eps == .005, K_hat],
          s[protocol == "rec" & rule == "adjacent_pointwise" & eps == .005, K_hat])
  exp_pw <- if (dn == "g100") c(50L, 80L, NA, NA) else c(50L, 80L, 120L, 120L)
  gate(identical(pw, exp_pw), "2d %s: pointwise selections reproduce Table S19 (%s)",
       dn, paste(pw, collapse = "/"))
}
OUT$mdna_sens <- sens
emit(sens[, .(design, c_part, n_models_in_support, n_coarse, protocol, rule, eps, K_hat,
              certified, n_comparisons, M, max_abs_dR2)], "mdna_selection_design_sensitivity")

# ================================== 2e ==========================================
f <- as.numeric(Matrix::colSums(prep$dtm_train))
f_str <- .freq_strata(prep$dtm_train, 5L)
rw <- copy(x$resid_words)[, stratum := f_str[match(word, colnames(prep$dtm_train))]]
mass <- rw[, .(words = .N, mass_pp = 100 * sum(resid_mean)), by = stratum][order(stratum)]
gate(near(mass$mass_pp, c(0.016842, 0.010720, 0.001031, -0.072319, 0.043726), 5e-7) &&
     abs(sum(mass$mass_pp)) < 1e-9 && near(sum(abs(mass$mass_pp)) / 2, 0.072319, 5e-7),
     "2e residual mass by training-frequency group (pp): %s; half abs sum %.6f",
     paste(sprintf("%+.6f", mass$mass_pp), collapse = " "), sum(abs(mass$mass_pp)) / 2)
OUT$mdna_mass_K50 <- mass
emit(mass, "mdna_residual_mass_K50_cache")

# MD&A comparators on a common footing: the protocol-matched likelihood optimum.
# On the harmonised cells and the retained documents, sum_j D_j(K) = 2 (C - l_K)
# with l_K the log score of the scored tokens, and the cells / retained set do
# not depend on K: argmax_K Micro-Deviance == argmin_K of that log loss.
opt_mdna <- x$summary[metric == "dev", .(K_opt_micro = K[which.max(r2_micro)],
                                        K_opt_macro = K[which.max(r2_macro)]), by = "eval"]
OUT$mdna_matched_opt <- opt_mdna
emit(opt_mdna, "mdna_matched_optimum")

# ================================== 2f ==========================================
sim_block <- function(file, label) {
  e1 <- qs_read(file); cfg <- e1$config; Kt <- cfg$K_true
  g <- e1$gains[metric == "dev"]
  # (i) simultaneous adjacent vs pointwise, from cached gains, all seeds
  sa <- rbindlist(lapply(split(g, by = c("eval", "replicate")), function(gg) {
    rbindlist(lapply(EPS, function(e) data.table(
      eval = gg$eval[1L], replicate = gg$replicate[1L], eps = e,
      pointwise = select_k_epsilon(gg, e, ALPHA)$K_hat,
      sim_adjacent = select_k_adjacent_simultaneous(gg, e, ALPHA)$K_hat)))
  }))
  # (ii) lower bound on the total-gain selection, all seeds: a certified K must
  # pass (a) its adjacent bound at the ALL-PAIRS quantile and (b) every point
  # estimate of the gain to a larger K <= eps (means of paired gains equal
  # differences of Macro means because the retained set is K-invariant).
  m <- length(cfg$K_grid); zM <- qnorm(1 - ALPHA / (m * (m - 1) / 2))
  mac <- e1$summary[metric == "dev" & eval %in% unique(g$eval),
                    .(eval, replicate, K, r2_macro)]
  lb <- rbindlist(lapply(split(g, by = c("eval", "replicate")), function(gg) {
    mm <- mac[eval == gg$eval[1L] & replicate == gg$replicate[1L]][order(K)]
    rbindlist(lapply(EPS, function(e) {
      ok <- vapply(seq_len(nrow(gg)), function(i) {
        k_i <- gg$K[i]
        # isTRUE: an undefined bound (collapsed support) certifies nothing
        isTRUE((gg$delta_mean[i] + zM * gg$se[i] <= e) &&
                 all(mm[K > k_i, r2_macro] - mm[K == k_i, r2_macro] <= e))
      }, NA)
      data.table(eval = gg$eval[1L], replicate = gg$replicate[1L], eps = e,
                 total_gain_lower_bound = if (any(ok)) min(gg$K[ok]) else NA_integer_)
    }))
  }))
  # (iii) exact rules for the one replicate with document scores
  ex <- rbindlist(lapply(c("ho_reconstruction", "ho_completion"), function(mn)
    select_k_all_rules(e1$doc_rep1[eval == mn & metric == "dev"], EPS, ALPHA)[
      , `:=`(eval = mn, replicate = 1L)]), use.names = TRUE)
  # (iv) protocol-matched likelihood optimum vs the reported perplexity comparator
  opt <- e1$summary[metric == "dev", .(K_opt = K[which.max(r2_micro)]),
                    by = c("eval", "replicate")]
  ppx <- e1$comparators[metric == "held_out_perplexity",
                        .(K_min = K[which.min(value)]), by = replicate]
  # The Micro-Deviance / log-score identity is exact on COMMON cells and
  # documents; the stored comparator is an UNBINNED perplexity over all
  # documents, so agreement with the reconstruction Micro optimum is expected
  # but not guaranteed seed by seed.
  agree <- merge(opt[eval == "ho_reconstruction", .(replicate, K_opt)], ppx,
                 by = "replicate")[, .(n = .N, seeds_agreeing = sum(K_opt == K_min))]
  tab <- merge(sa, lb, by = c("eval", "replicate", "eps"))
  summ <- tab[, .(n = .N,
    pointwise_exact = mean(hit(pointwise, Kt)), pointwise_mean = mean(pointwise, na.rm = TRUE),
    pointwise_na = sum(is.na(pointwise)),
    sim_adj_exact = mean(hit(sim_adjacent, Kt)), sim_adj_mean = mean(sim_adjacent, na.rm = TRUE),
    sim_adj_na = sum(is.na(sim_adjacent)),
    tg_lb_exact = mean(hit(total_gain_lower_bound, Kt)),
    tg_lb_mean = mean(total_gain_lower_bound, na.rm = TRUE),
    tg_lb_na = sum(is.na(total_gain_lower_bound))), by = c("eval", "eps")]
  optd <- opt[, .(n = .N, exact = mean(K_opt == Kt), mean_K = mean(K_opt),
                  dist = paste(sprintf("%d:%d", sort(unique(K_opt)),
                                       as.integer(table(K_opt))), collapse = " ")),
              by = "eval"]
  # (v) the other discrepancy families under both adjacent rules (cached gains)
  fam <- rbindlist(lapply(split(e1$gains, by = c("metric", "eval", "replicate")),
    function(gg) rbindlist(lapply(EPS, function(e) data.table(
      metric = gg$metric[1L], eval = gg$eval[1L], replicate = gg$replicate[1L],
      eps = e, pointwise = select_k_epsilon(gg, e, ALPHA)$K_hat,
      sim_adjacent = select_k_adjacent_simultaneous(gg, e, ALPHA)$K_hat)))))
  mode_of <- function(v) { v <- v[!is.na(v)]
    if (length(v)) as.integer(names(which.max(table(v)))) else NA_integer_ }
  fam_s <- fam[, .(n = .N, pw_mode = mode_of(pointwise), pw_exact = mean(hit(pointwise, Kt)),
                   pw_mean = mean(pointwise, na.rm = TRUE),
                   sa_mode = mode_of(sim_adjacent), sa_exact = mean(hit(sim_adjacent, Kt)),
                   sa_mean = mean(sim_adjacent, na.rm = TRUE),
                   sa_na = sum(is.na(sim_adjacent))), by = c("metric", "eval", "eps")]
  # (vi) share of documents below the null-discrepancy floor, by protocol
  excl <- e1$summary[K == Kt, .(excl_mean = mean(null_excl_share),
                                excl_min = min(null_excl_share),
                                excl_max = max(null_excl_share)),
                     by = c("metric", "eval")]
  list(label = label, K_true = Kt, per_seed = tab, summary = summ, exact_rep1 = ex,
       families = fam_s, exclusion = excl,
       matched_opt = opt, matched_opt_dist = optd,
       perplexity_dist = ppx[, .N, by = K_min][order(K_min)],
       perplexity_vs_micro_rec = agree)
}
E1 <- sim_block(p_data("E1", "e1_results_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"), "baseline")
s1 <- E1$summary
gate(near(s1[eval == "ho_completion" & eps == .01, sim_adj_exact], 1) &&
     near(s1[eval == "ho_reconstruction" & eps == .01, sim_adj_exact], 0.9) &&
     near(s1[eval == "ho_completion" & eps == .005, sim_adj_exact], 0.8),
     "2f E1 simultaneous adjacent: K* in %.0f/10 (com), %.0f/10 (rec) at .01; %.0f/10 (com) at .005",
     10 * s1[eval == "ho_completion" & eps == .01, sim_adj_exact],
     10 * s1[eval == "ho_reconstruction" & eps == .01, sim_adj_exact],
     10 * s1[eval == "ho_completion" & eps == .005, sim_adj_exact])
gate(E1$matched_opt_dist[eval == "ho_completion", dist] == "40:3 50:4 60:1 70:1 90:1" &&
     E1$matched_opt_dist[eval == "ho_reconstruction", dist] == "90:1 100:9" &&
     identical(E1$perplexity_dist$N, c(1L, 9L)),
     "2f matched optimum: completion [%s]; reconstruction [%s], same DISTRIBUTION as the reported perplexity comparator [90:1 100:9]",
     E1$matched_opt_dist[eval == "ho_completion", dist],
     E1$matched_opt_dist[eval == "ho_reconstruction", dist])
D2 <- sim_block(p_data("E1", "e1_results_E1_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01.qs2"), "second DGP")
OUT$sim <- list(baseline = E1, dgp2 = D2)
for (o in list(E1, D2)) {
  tg <- gsub(" ", "_", o$label)
  emit(o$summary, sprintf("sim_selection_rules_%s", tg))
  emit(o$per_seed, sprintf("sim_selection_rules_by_seed_%s", tg))
  emit(o$exact_rep1[, .(eval, rule, eps, K_hat, certified, n_comparisons, max_ub, M)],
       sprintf("sim_selection_exact_replicate1_%s", tg))
  emit(o$matched_opt_dist, sprintf("sim_matched_optimum_%s", tg))
  emit(o$families, sprintf("sim_selection_rules_by_family_%s", tg))
  emit(o$exclusion, sprintf("sim_floor_exclusion_%s", tg))
}

# ================================== 2g ==========================================
# The "truth" in the coverage experiment is the mean of ONE independent
# reference sample of J_truth documents per training fit, not a known number;
# its standard error was not saved. Recover the per-document SD from the
# evaluation SEs (reconstruction retains every document, so n = J_ev) and form
# (a) the coverage a PERFECTLY calibrated interval would show against a target
#     this noisy, 2 Phi(z / sqrt(1 + J_ev / J_truth)) - 1, and
# (b) a two-sample calibration check that adds the reference variance.
# (b) is a joint check of estimate and reference, NOT coverage of a known
# parameter. Monte Carlo error is computed across TRAINING FITS, because the
# ten evaluation draws of a fit share that fit and its reference sample.
cal_block <- function(file, label) {
  e2 <- qs_read(file); Jt <- e2$config$J_truth
  one <- function(dt, est, tru) {
    d <- copy(dt)
    d[, sd_unit := se * sqrt(J_ev)]
    d[, se_ref := mean(sd_unit) / sqrt(Jt), by = c("train_seed", "K")]
    d[, `:=`(cov_raw = lwr <= get(tru) & get(tru) <= upr,
             cov_joint = abs(get(est) - get(tru)) <= qnorm(.975) * sqrt(se^2 + se_ref^2))]
    fit <- d[, .(raw = mean(cov_raw), joint = mean(cov_joint)), by = c("J_ev", "train_seed")]
    fit[, .(coverage_raw = mean(raw), mc_se_raw = sd(raw) / sqrt(.N),
            coverage_joint = mean(joint), mc_se_joint = sd(joint) / sqrt(.N),
            min_fit = min(raw), max_fit = max(raw), n_fits = .N), by = J_ev][
      , expected_raw_if_calibrated := 2 * pnorm(qnorm(.975) / sqrt(1 + J_ev / Jt)) - 1][
      , se_ref_over_se_ev := sqrt(J_ev / Jt)][order(J_ev)]
  }
  fitc <- one(e2$cover, "r2_macro", "mu")[, target := "average fit (Macro)"]
  gapc <- one(e2$gap, "gap", "gap_true")[, target := "Micro-Macro gap"]
  # selection by evaluation size under the simultaneous adjacent rule
  g <- copy(e2$gains)[, metric := "dev"]
  selJ <- rbindlist(lapply(split(g, by = c("train_seed", "J_ev", "rep")), function(gg)
    rbindlist(lapply(EPS, function(e) data.table(
      J_ev = gg$J_ev[1L], eps = e,
      pointwise = { i <- which(!is.na(gg$ub_onesided) & gg$ub_onesided <= e)
                    if (length(i)) min(gg$K[i]) else NA_integer_ },
      sim_adjacent = select_k_adjacent_simultaneous(gg, e, ALPHA)$K_hat)))))
  Kt <- e2$config$K_true
  selS <- selJ[, .(n = .N, pointwise_exact = mean(hit(pointwise, Kt)),
                   pointwise_na = sum(is.na(pointwise)),
                   sim_adj_exact = mean(hit(sim_adjacent, Kt)),
                   sim_adj_mean = mean(sim_adjacent, na.rm = TRUE),
                   sim_adj_na = sum(is.na(sim_adjacent))), by = c("J_ev", "eps")][order(eps, J_ev)]
  list(label = label, J_truth = Jt, cover_tab = e2$cover_tab,
       calibration = rbindlist(list(fitc, gapc)), selection_by_Jev = selS)
}
C1 <- cal_block(p_data("E2", "e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"), "baseline")
cf <- C1$calibration[target == "average fit (Macro)"]
gate(near(cf$coverage_raw, c(0.960, 0.950, 0.945), 5e-4) &&
     near(range(C1$cover_tab$coverage), c(0.90, 0.99), 5e-3) && C1$J_truth == 2000,
     "2g raw fit coverage by J_ev = %s (cells %.2f-%.2f, as in Table S3); J_truth = %d",
     paste(sprintf("%.3f", cf$coverage_raw), collapse = "/"),
     min(C1$cover_tab$coverage), max(C1$cover_tab$coverage), C1$J_truth)
gate(near(C1$selection_by_Jev[eps == .01, pointwise_exact], c(0.57, 0.76, 0.88), 5e-3),
     "2g pointwise Pr(K-hat = K*) by J_ev reproduces the paper: %s",
     paste(sprintf("%.2f", C1$selection_by_Jev[eps == .01, pointwise_exact]), collapse = "/"))
C2 <- cal_block(p_data("E2", "e2_results_E2_full_Kstar20_J1000_W10000_a0p5_b0p01_k5-50by5_warplda_fa0p1_fb0p01.qs2"), "second DGP")
OUT$cal <- list(baseline = C1, dgp2 = C2)
for (o in list(C1, C2)) {
  tg <- gsub(" ", "_", o$label)
  emit(o$calibration, sprintf("sim_coverage_calibration_%s", tg))
  emit(o$selection_by_Jev, sprintf("sim_selection_by_Jev_%s", tg))
}

# ================================== 2h ==========================================
# Planted-word recovery (Study III) under the word-level null deviance the paper
# DEFINES. OpTop <= 0.20.1 returns the word-level null deviance without its
# Poisson linear term, 2 sum N log(N/B), while the fitted deviance carries it.
# Held-out (training baseline, evaluation counts) the omitted term 2 (N_w - B_w)
# is not zero, so held-out word-level R2 and every ranking built on it change.
# d_model and B_w are cached per word; only N_w is needed, from the evaluation
# corpus regenerated deterministically from the stored seeds (no fit, no fold-in).
source(here::here("Code", "config", "configs.R"))
f_e4 <- p_data("E4", "e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2")
e4 <- qs_read(f_e4)
c4 <- e4$config; sd4 <- make_seeds(c4$seed_base, 1L)
tr4 <- sim_lda_corpus(c4$J_train, c4$W, c4$K_true, c4$alpha_DGP, c4$beta_DGP,
                      c4$length_spec, seed = sd4$dgp_seed, doc_prefix = "tr")
s_max <- max(c4$alternatives$contamination_eval$strengths)
pi_stop4 <- make_stopword_dist(c4$W, c4$n_stop, seed = sd4$dgp_seed + 901L)
ev4 <- sim_lda_corpus(c4$J_ev, c4$W, c4$K_true, c4$alpha_DGP, c4$beta_DGP,
                      c4$length_spec, seed = sd4$eval_seed + 4000017L, Phi = tr4$Phi,
                      contamination = list(pi_stop = pi_stop4, w_mean = s_max,
                                           w_conc = 10, mode = "shared"),
                      doc_prefix = "ev")
wr <- copy(e4$word_rank)[match(colnames(ev4$dtm), word_id)]
N_w4 <- as.numeric(Matrix::colSums(ev4$dtm))
pi4 <- as.numeric(OpTop::optop_make_baseline(tr4$dtm)$pi_glob)
gate(near(wr$B_w, pi4 * sum(ev4$dtm), 1e-8) &&
     identical(wr$planted, pi_stop4 > 0) &&
     near(wr$doc_freq, as.numeric(Matrix::colSums(ev4$dtm > 0)), 0),
     "2h regenerated evaluation corpus reproduces the cached word table (B_w, doc_freq, planted flags)")
# The linear term is subtracted ONLY from a word table written under the legacy
# convention. The convention is read from the object's scoring stamp or, for the
# frozen pre-revision cache, established by its content hash; anything else stops
# here (word_null_convention_of). A table already in Poisson form passes through.
conv4 <- word_null_convention_of(e4, f_e4)
log_msg("2h word-null convention of the E4 input: %s", conv4)
corrected4 <- revision_word_null(wr, N_w4, wr$B_w, conv4)
wr[, d_null_paper := corrected4$d_null]
wr[, r2_paper := ifelse(d_null_paper > 0, 1 - d_model / d_null_paper, NA_real_)]
rank_tab <- function(col, label) {
  k <- wr[keep == TRUE & !is.na(get(col))][order(get(col))]
  np <- sum(k$planted)
  data.table(null_deviance = label, n_words = nrow(k), n_planted = np,
             precision_at_n = mean(head(k$planted, np)),
             recall_top50 = sum(head(k$planted, 50L)) / np,
             median_rank_planted = median(which(k$planted)),
             median_r2_planted = median(k[[col]][k$planted]),
             median_r2_other = median(k[[col]][!k$planted]))
}
pl <- rbindlist(list(rank_tab("r2_word", "package: no linear term"),
                     rank_tab("r2_paper", "paper: Poisson form")))
gate(pl[1L, n_planted] == 37L && near(pl[1L, median_rank_planted], 2716, 0.5),
     "2h package convention reproduces Table S10 (37 planted words kept, median rank %.0f)",
     pl[1L, median_rank_planted])
OUT$e4_planted <- pl
emit(pl, "e4_planted_words_null_conventions")

f_out <- p_data("MDNA", sprintf("revision_postprocess%s.qs2", SFX))
cache_put(OUT, f_out, list(experiment = "REVISION", profile = "postprocess"))
log_msg("saved %s", basename(f_out))

cat("\n================ SUMMARY ================\n")
cat("\nMD&A Deviance selections (document-level SE):\n")
print(dcast(sel[metric == "dev"], rule ~ protocol + eps, value.var = "K_hat"))
cat("\nMD&A, all families, completion:\n")
print(dcast(sel[protocol == "com"], rule ~ metric + eps, value.var = "K_hat"))
cat("\nFirms absent from training (firm-cluster SE):\n")
print(dcast(sel_abs[metric == "dev"], rule ~ protocol + eps, value.var = "K_hat"))
cat("\ndelta sensitivity (c = 1):\n"); print(dl_w)
cat("\nDesign sensitivity on the refined-support caches:\n")
print(dcast(sens, design + n_models_in_support + rule ~ protocol + eps, value.var = "K_hat"))
cat("\nMD&A matched optimum:\n"); print(opt_mdna)
for (o in list(E1, D2)) {
  cat(sprintf("\n--- simulation: %s (K* = %d) ---\n", o$label, o$K_true))
  print(o$summary); print(o$matched_opt_dist)
  cat("families (adjacent rules):\n"); print(o$families[eps == 0.01])
  cat("floor exclusion share at K*:\n"); print(o$exclusion)
  cat("reported perplexity comparator:\n"); print(o$perplexity_dist)
  cat("seeds where it shares the reconstruction Micro-Deviance optimum:\n")
  print(o$perplexity_vs_micro_rec)
  cat("exact rules, replicate 1:\n")
  print(dcast(o$exact_rep1, rule ~ eval + eps, value.var = "K_hat"))
}
cat("\nStudy III planted words, word-level fit ranking under the two null conventions:\n")
print(pl, digits = 4)
for (o in list(C1, C2)) {
  cat(sprintf("\n--- coverage calibration: %s (J_truth = %d) ---\n", o$label, o$J_truth))
  print(o$calibration, digits = 3); print(o$selection_by_Jev, digits = 3)
}
log_msg("=== postprocess_revision complete (suffix %s) ===", SFX)
