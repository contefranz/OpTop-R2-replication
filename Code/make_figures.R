# =============================================================================
# make_figures.R -- all Section 5 figures (F1-F12) from saved results only.
# Usage:  Rscript Code/make_figures.R [smoke|pilot|full] [exp=E3,E6] [label=name]
#   exp=   restrict to a subset of experiments (default: all)
#   label= output tree Results/Figures/<label>/ (default: profile[+overrides])
# Skips any figure whose inputs are missing; never recomputes.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))
cli <- parse_cli(args)
PROFILE <- cli$profile
# experiment selector: strip exp= before it can reach the label or apply_overrides
sel_exp <- setdiff(toupper(strsplit(cli$overrides$exp %||% "", ",")[[1]]), "")
cli$overrides$exp <- NULL
RUN_LABEL <- run_label(list(profile = PROFILE), cli$overrides, cli$label)
log_msg("figures for run label '%s'%s", RUN_LABEL,
        if (length(sel_exp)) sprintf(" (experiments: %s)",
                                     paste(sel_exp, collapse = ",")) else "")

load_results <- function(exp) {
  if (length(sel_exp) && !toupper(exp) %in% sel_exp) return(NULL)
  cfg <- apply_overrides(get_config(exp, PROFILE), cli$overrides, strict = FALSE)
  f <- p_data(exp, sprintf("%s_results_%s.qs2", tolower(exp), run_tag(cfg)))
  if (!file.exists(f)) { log_msg("skip %s (no results: %s)", exp, basename(f)); return(NULL) }
  cache_get(f)
}

e1 <- load_results("E1"); e2 <- load_results("E2"); e3 <- load_results("E3")
e4 <- load_results("E4"); e5 <- load_results("E5"); e6 <- load_results("E6")

# ============ F1: held-out fit curves, all three discrepancy families =============
if (!is.null(e1)) {
  K_true <- e1$config$K_true
  add_labs <- function(dt) dt[, `:=`(
    eval_lab = factor(eval_label(eval), levels = EVAL_LEVELS),
    metric_lab = factor(metric_label(metric), levels = METRIC_LEVELS))]
  s <- add_labs(copy(e1$summary))
  band <- s[, .(lo = min(r2_macro), hi = max(r2_macro),
                lo_mi = min(r2_micro), hi_mi = max(r2_micro)),
            by = .(K, eval_lab, metric_lab)]
  ci1 <- add_labs(copy(e1$ci)[replicate == 1L])
  s1 <- s[replicate == 1L]

  p1 <- ggplot(s1, aes(x = K, y = r2_macro, linetype = eval_lab)) +
    geom_ribbon(data = band, aes(x = K, ymin = lo, ymax = hi, group = eval_lab),
                inherit.aes = FALSE, fill = "grey85", alpha = 0.5) +
    geom_ribbon(data = ci1, aes(x = K, ymin = lwr, ymax = upr, group = eval_lab),
                inherit.aes = FALSE, fill = "grey60", alpha = 0.35) +
    geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.7) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    scale_eval_linetype() + scale_eval_shape() +
    facet_wrap(~ metric_lab, scales = "free_y") +
    labs(x = "Number of topics (K)", y = expression(R[Macro]^2),
         title = "Held-out fit over K by discrepancy family",
         subtitle = "Deviance is the primary index; Pearson and Squared-Error shown for comparison. Ribbons: Prop.-2 95% CI (one seed) + across-seed range.") +
    theme_paper()
  save_fig(p1, "F1_heldout_curves", "E1", width = 11, height = 4.6)

  p1b <- ggplot(s1, aes(x = K, y = r2_micro, linetype = eval_lab)) +
    geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.7) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    scale_eval_linetype() + scale_eval_shape() +
    facet_wrap(~ metric_lab, scales = "free_y") +
    labs(x = "Number of topics (K)", y = expression(R[Micro]^2),
         title = "Micro aggregation by discrepancy family (companion to F1)") +
    theme_paper()
  save_fig(p1b, "F1b_micro_curves", "E1", width = 11, height = 4.6)
}

# ============= F2: adjacent gains and the epsilon rule, all families ==============
if (!is.null(e1)) {
  gm_labs <- function(dt) dt[, `:=`(
    eval_lab = droplevels(factor(eval_label(eval), levels = EVAL_LEVELS)),
    metric_lab = factor(metric_label(metric), levels = METRIC_LEVELS))]
  g1 <- gm_labs(copy(e1$gains)[replicate == 1L])
  kh <- gm_labs(copy(e1$khat)[replicate == 1L & rule == "eps_0p01" &
                                metric %in% c("dev", "chisq", "se")])
  p2 <- ggplot(g1, aes(x = K, y = delta_mean)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_hline(yintercept = c(0.01, 0.005), linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = lwr, ymax = ub_onesided), width = 0.3, colour = "grey25") +
    geom_point(size = 1.5) +
    geom_vline(data = kh[!is.na(K_hat)], aes(xintercept = K_hat),
               linetype = "dotdash", colour = "grey20") +
    facet_grid(eval_lab ~ metric_lab, scales = "free_y") +
    labs(x = "Number of topics (K); gain to the next grid point",
         y = expression(Delta * R^2),
         title = "Held-out adjacent gains with one-sided upper bounds",
         subtitle = "Deviance is the selection index; dashed: eps = 0.01, 0.005; dot-dash: selected K-hat (eps = 0.01)") +
    theme_paper()
  save_fig(p2, "F2_adjacent_gains", "E1", width = 11, height = 6)
}

# =================== F3: CI coverage and t-statistic QQ (E2) =======================
if (!is.null(e2)) {
  ct <- e2$cover_tab
  p3a <- ggplot(ct, aes(x = K, y = coverage)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40") +
    geom_line(colour = "grey20") + geom_point(size = 1.6) +
    facet_wrap(~J_ev, labeller = label_both) +
    coord_cartesian(ylim = c(0.8, 1)) +
    labs(x = "Number of topics (K)", y = "Empirical coverage",
         title = "Coverage of the 95% held-out Macro CI (Prop. 2)") +
    theme_paper()
  ts <- e2$cover[K == e2$config$K_true]
  p3b <- ggplot(ts, aes(sample = t_stat)) +
    stat_qq(size = 0.7, alpha = 0.5) + stat_qq_line(colour = "grey40") +
    facet_wrap(~J_ev, labeller = label_both) +
    labs(x = "N(0,1) quantiles", y = "t-statistic quantiles",
         title = sprintf("t-statistics at K = %d", e2$config$K_true)) +
    theme_paper()
  save_fig(p3a / p3b, "F3_coverage_qq", "E2", width = 8.6, height = 8.4)
}

# ==================== F4: moment-test power curves (E4) ============================
if (!is.null(e4) && nrow(e4$power)) {
  K_true <- e4$config$K_true
  pw <- e4$power_tab[K == K_true]
  sz <- e4$size_tab[K == K_true]
  if (!"size_raw" %in% names(sz)) sz[, size_raw := size]  # pre-revision objects
  # a line needs >=2 strengths; single-strength alternatives render as points
  # only (.SD[cond] keeps the columns even when EVERY alt is single-strength,
  # e.g. the smoke profile -- `if (cond) .SD` would drop them and break aes())
  pw_line <- pw[, .SD[uniqueN(strength) > 1], by = alt]
  p4 <- ggplot(pw, aes(x = strength, y = power, colour = test, group = test)) +
    geom_hline(data = sz, aes(yintercept = size_raw, colour = test),
               linetype = "dotted") +
    geom_hline(yintercept = 0.05, colour = "grey70") +
    geom_line(data = pw_line, linewidth = 0.8) + geom_point(size = 1.8) +
    facet_wrap(~alt, scales = "free_x") +
    scale_colour_grey(start = 0, end = 0.6, name = "Test") +
    labs(x = "Misspecification strength", y = "Rejection rate (5% level)",
         title = sprintf("Power of the moment tests at K = %d", K_true),
         subtitle = "Dotted: raw rejection under the correct DGP (conditional null); grey: nominal 5%") +
    theme_paper()
  save_fig(p4, "F4_power_curves", "E4", width = 9.2)
}

# ============ F5: "fit barely moves, tests fire" + planted-word ranks ==============
if (!is.null(e4) && nrow(e4$r2)) {
  K_true <- e4$config$K_true
  r2a <- e4$r2_tab[K == K_true]
  pw <- e4$power_tab[K == K_true, .(power = mean(power)), by = .(alt, strength)]
  both <- merge(r2a, pw, by = c("alt", "strength"))
  b_long <- melt(both, id.vars = c("alt", "strength"),
                 measure.vars = c("r2_macro", "power"),
                 variable.name = "what")
  b_long[, what := fifelse(what == "r2_macro",
                           "Held-out Macro R2 (Dev)", "Moment-test power")]
  bl_line <- b_long[, .SD[uniqueN(strength) > 1], by = alt]
  p5a <- ggplot(b_long, aes(x = strength, y = value, linetype = what)) +
    geom_line(data = bl_line, linewidth = 0.8) + geom_point(size = 1.8) +
    facet_wrap(~alt, scales = "free_x") +
    labs(x = "Misspecification strength", y = NULL, linetype = NULL,
         title = "Overall fit vs specification tests under misspecification") +
    theme_paper()
  save_fig(p5a, "F5_fit_vs_tests", "E4", width = 9.2)

  if (!is.null(e4$word_resid)) {
    wres <- e4$word_resid[order(-resid_mean)]
    wres[, rank := .I]
    p5b <- ggplot(wres[rank <= 300], aes(x = rank, y = resid_mean,
                                         colour = planted)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_point(size = 1.1, alpha = 0.8) +
      scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "black"),
                          labels = c("regular", "planted"), name = NULL) +
      labs(x = "Rank (largest mean held-out residual first)",
           y = expression(bar(e)[w] ~ "(probability mass)"),
           title = "Word-level residual diagnostic: planted words surface",
           subtitle = "Eval-only contamination; mean residual = observed - fitted word probability") +
      theme_paper()
    save_fig(p5b, "F5b_word_ranks", "E4")
  }
}

# ================= F6: Micro-Macro gap and its decomposition (E3) ==================
if (!is.null(e3)) {
  K_true <- e3$config$K_true
  sc_lab <- as_labeller(sapply(e3$config$scenarios, `[[`, "label"))
  gp <- copy(e3$gap)[eval == "ho_reconstruction"]
  gp[, metric_lab := factor(metric_label(metric), levels = METRIC_LEVELS)]
  gp1 <- gp[replicate == 1L]
  p6a <- ggplot(gp1, aes(x = K, y = gap)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey80", alpha = 0.6) +
    geom_line(data = gp, aes(group = replicate), colour = "grey55",
              linewidth = 0.3) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    facet_grid(metric_lab ~ scenario, scales = "free_y",
               labeller = labeller(scenario = sc_lab, metric_lab = label_value)) +
    labs(x = "Number of topics (K)",
         y = expression(R["Micro"]^2 - R["Macro"]^2),
         title = "Held-out Micro-Macro gap by discrepancy family, with delta-method 95% CI",
         subtitle = "Deviance primary; thin lines: other replicates") +
    theme_paper()
  save_fig(p6a, "F6_gap_ci", "E3", width = 10.5, height = 7)

  dc <- e3$decomp[eval == "ho_reconstruction" & replicate == 1L]
  dcl <- melt(dc, id.vars = c("K", "scenario"),
              measure.vars = c("ch_length", "ch_atypicality", "ch_interaction"),
              variable.name = "channel")
  dcl[, channel := factor(channel,
                          levels = c("ch_length", "ch_atypicality", "ch_interaction"),
                          labels = c("Length", "Atypicality", "Interaction"))]
  p6b <- ggplot(dcl, aes(x = K, y = value, fill = channel)) +
    geom_col(position = "stack", width = 0.8) +
    geom_line(data = dc, aes(x = K, y = gap), inherit.aes = FALSE,
              linewidth = 0.7) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    scale_fill_grey(start = 0.2, end = 0.85, name = "Channel") +
    facet_wrap(~scenario, labeller =
                 as_labeller(sapply(e3$config$scenarios, `[[`, "label"))) +
    labs(x = "Number of topics (K)", y = "Gap contribution",
         title = "Exact channel decomposition of the gap (Corollary 1)",
         subtitle = "Solid line: total gap") +
    theme_paper()
  save_fig(p6b, "F6b_gap_decomposition", "E3", width = 9.6)
}

# =================== F7: document-level heterogeneity scatters =====================
if (!is.null(e3) && nrow(e3$scatter)) {
  scA <- e3$scatter[scenario == "A"]
  scC <- e3$scatter[scenario == "C"]
  p7a <- ggplot(scA, aes(x = L, y = r2_doc)) +
    geom_point(aes(shape = group), alpha = 0.35, size = 1.1) +
    stat_summary_bin(fun = median, bins = 12, geom = "line",
                     colour = "black", linewidth = 0.8) +
    scale_x_log10() +
    labs(x = "Document length (tokens, log scale)",
         y = expression(R["Dev,j"]^2),
         title = "A: fit vs length", shape = NULL) +
    theme_paper()
  p7b <- ggplot(scC, aes(x = kappa, y = r2_doc)) +
    geom_point(aes(shape = group), alpha = 0.35, size = 1.1) +
    stat_summary_bin(fun = median, bins = 12, geom = "line",
                     colour = "black", linewidth = 0.8) +
    labs(x = expression(kappa[j] ~ "(KL divergence from baseline)"),
         y = expression(R["Dev,j"]^2),
         title = "C: fit vs atypicality", shape = NULL) +
    theme_paper()
  save_fig(p7a + p7b, "F7_doc_scatters", "E3", width = 10.2, height = 4.6)
}

# F8 (estimator robustness: Gibbs / multi-start dispersion) was removed in the
# 4-study consolidation together with E1's robustness arms.

# ========================= F9: design sensitivity (E5) =============================
if (!is.null(e5)) {
  K_true <- e5$config$K_true
  p9a <- ggplot(e5$grid_curves[metric == "dev"],
                aes(x = K, y = r2_micro, colour = grid)) +
    geom_line(linewidth = 0.8) + geom_point(aes(shape = grid), size = 1.7) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    scale_colour_grey(start = 0, end = 0.65, name = "Estimation grid") +
    scale_shape(name = "Estimation grid") +
    labs(x = "Number of topics (K)", y = expression(R["Dev,Micro"]^2),
         title = "Grid-extension sensitivity of the harmonized support") +
    theme_paper()
  cc <- e5$c_curves
  p9b <- ggplot(cc, aes(x = K, y = r2_micro,
                        linetype = factor(c_value))) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = factor(c_value)), size = 1.7) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    facet_wrap(~metric, labeller = as_labeller(metric_label)) +
    scale_shape(name = "c") +
    labs(x = "Number of topics (K)", y = expression(R^2),
         linetype = "c",
         title = "Threshold sensitivity: c = 1 vs c = 5") +
    theme_paper()
  save_fig(p9a + p9b, "F9_design_sensitivity", "E5", width = 11, height = 4.6)
}

# ========= F10 (appendix): E1 robustness variants (K* = 20; alpha = 0.1) ==========
for (vex in c("E1b", "E1c")) {
  ev <- load_results(vex)
  if (is.null(ev)) next
  sv <- ev$summary[metric == "dev" & eval != "insample"]
  sv[, eval_lab := factor(eval_label(eval), levels = EVAL_LEVELS)]
  sv_m <- sv[, .(r2_macro = mean(r2_macro), lo = min(r2_macro),
                 hi = max(r2_macro)), by = .(K, eval_lab)]
  pv <- ggplot(sv_m, aes(x = K, y = r2_macro, linetype = eval_lab)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = eval_lab),
                fill = "grey85", alpha = 0.5, colour = NA) +
    geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.7) +
    geom_vline(xintercept = ev$config$K_true, linetype = "dashed",
               colour = "grey40") +
    scale_eval_linetype() + scale_eval_shape() +
    labs(x = "Number of topics (K)", y = expression(R["Dev,Macro"]^2),
         title = sprintf("Robustness variant %s: K* = %d, alpha = %s",
                         vex, ev$config$K_true,
                         # DGP concentration; fall back for pre-rename caches
                         ev$config$alpha_DGP %||% ev$config$alpha),
         subtitle = "Band: across-seed range of the seed means") +
    theme_paper()
  save_fig(pv, sprintf("F10_%s_curves", vex), "E1")
}

# ============= F11/F12: word-level dual perspective (E6, §3.8) =====================
if (!is.null(e6)) {
  K_true <- e6$config$K_true
  # Paper uses the correctly specified scenario only: the former W2
  # (shared both-sides stopwords) is absorbed by refitting and was dropped
  # from the study; its data remains in the results object.

  # F11: w-Micro and w-Macro over K, held-out (W1); gap shaded
  cv <- e6$curve[eval == "ho_reconstruction" & scenario == "W1",
                 .(w_Micro = mean(r2_micro_word), w_Macro = mean(r2_macro_word)),
                 by = K]
  cl <- melt(cv, id.vars = "K",
             measure.vars = c("w_Micro", "w_Macro"),
             variable.name = "agg", value.name = "r2")
  p11 <- ggplot(cl, aes(x = K, y = r2, linetype = agg)) +
    geom_ribbon(data = cv, aes(x = K, ymin = w_Macro, ymax = w_Micro),
                inherit.aes = FALSE, fill = "grey85", alpha = 0.6) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    geom_vline(xintercept = K_true, linetype = "dashed", colour = "grey40") +
    scale_linetype_manual(values = c(w_Micro = "solid", w_Macro = "dotted"),
                          name = "Aggregation") +
    labs(x = "Number of topics (K)", y = expression(R["Dev,w"]^2),
         title = "Word-level dual perspective: w-Micro vs w-Macro over K",
         subtitle = "Shaded = w-Micro - w-Macro gap (frequency-space twin of the length-bias gap)") +
    theme_paper()
  save_fig(p11, "F11_word_micro_macro", "E6", width = 7.4)

  # F12: per-word held-out R^2 at K* vs document frequency (W1)
  if (nrow(e6$word_star)) {
    ws <- e6$word_star[keep == TRUE & scenario == "W1"]
    if (nrow(ws)) {
      p12 <- ggplot(ws, aes(x = log10_docfreq, y = r2_word)) +
        geom_point(size = 0.9, alpha = 0.5, colour = "grey55") +
        geom_smooth(se = FALSE, colour = "black", linewidth = 0.7,
                    method = "loess", formula = y ~ x) +
        labs(x = expression(log[10]("document frequency + 1")),
             y = bquote(R["Dev,w"]^2 ~ "at K =" ~ .(K_true)),
             title = "Word-level held-out fit vs frequency",
             subtitle = "Fit concentrates in the common vocabulary; the rare tail is estimation-hard") +
        theme_paper()
      save_fig(p12, "F12_word_fit_vs_freq", "E6", width = 7.4)
    }
  }
}

log_msg("figures complete (profile %s)", PROFILE)
