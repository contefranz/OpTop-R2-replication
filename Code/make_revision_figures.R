# =============================================================================
# make_revision_figures.R -- the figures whose CONTENT (not just caption) the
# Sept-2026 revision changes. Built from saved results only; written to NEW
# trees Results/Figures/<label>_rev1/, so no existing PDF is overwritten.
#
#   MD&A  R1  fit curves, the three completion selections marked       (main Fig. 2)
#   MD&A  R4  moment diagnostics with K-specific Test-3 strata, at the
#             reference fit and at the total-gain selection            (main Fig. 3)
#   E1    F2  adjacent gains; old "selected K-hat" marker replaced by the
#             simultaneous adjacent and total-gain selections          (supp Fig. S1)
#   E4    F4  raw rejection curves, Tests 1-2 only; subtitle no longer calls
#             baseline rejection "size"                                (supp Fig. S8)
#   E4    F5  fit vs rejection, Tests 1-2 only                         (supp Fig. S9)
# Test 3 is omitted from F4/F5 because its K* = 40 values were computed with the
# wrong strata and await the deferred E4 batch.
#
# Usage:  Rscript Code/make_revision_figures.R [suffix=_rev1]
# =============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
.args <- commandArgs(trailingOnly = TRUE)
SFX <- sub("^suffix=", "", grep("^suffix=", .args, value = TRUE))
if (!length(SFX)) SFX <- "_rev1"
L_MDNA <- paste0("mdna_2015_2016", SFX); L_SIM <- paste0("warptestFULL_NULL", SFX)
ALPHA <- 0.05

# ------------------------------- MD&A -------------------------------------------
x <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016.qs2"))
pp <- qs_read(p_data("MDNA", sprintf("revision_postprocess%s.qs2", SFX)))
rs <- qs_read(p_data("MDNA", sprintf("mdna_rescore%s.qs2", SFX)))
sel <- pp$mdna_sel[metric == "dev" & protocol == "com" & eps == 0.01]
marks <- data.table(
  K = c(sel[rule == "adjacent_pointwise", K_hat], sel[rule == "adjacent_simultaneous", K_hat],
        sel[rule == "total_gain", K_hat]),
  lab = c("pointwise", "simult. adjacent", "total gain"),
  lt = c("dotted", "dashed", "dotdash"))
eval_lab_m <- c(ins = "In-sample", rec = "Held-out (reconstruction)",
                com = "Held-out (completion)")
s <- copy(x$summary)[metric == "dev"]
s[, eval_lab := factor(eval_lab_m[eval], levels = EVAL_LEVELS)]
y_bot <- min(s$r2_macro)     # the lower part of the panel is free of curves at
                             # all three selections, so the labels sit there
p1 <- ggplot(s, aes(x = K, y = r2_macro, linetype = eval_lab)) +
  geom_vline(data = marks, aes(xintercept = K), colour = "grey55", linewidth = 0.5) +
  geom_text(data = marks, aes(x = K, y = y_bot, label = sprintf("%s: %d", lab, K)),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = -0.5, size = 3.1,
            colour = "grey20") +
  geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.5) +
  scale_eval_linetype() + scale_eval_shape() +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Number of topics (K)", y = expression(R["Macro"]^2),
       title = "MD&A corpus (fiscal 2015-2016): Deviance fit over K",
       subtitle = "Vertical lines: completion selections at eps = 0.01 under the three rules") +
  theme_paper()
save_fig(p1, "R1_fit_curves", "MDNA", width = 7.4, height = 4.6, label = L_MDNA)

K_ref <- as.integer(rs$params$K_ref); K_sel <- as.integer(rs$params$K_sel)
tl <- c(T1_freq_contrast = "Test 1: frequency contrast", T2_freq_strata = "Test 2: frequency strata",
        T3_fit_strata = "Test 3: training-fit strata")
st <- rs$tests$strata[K %in% c(K_ref, K_sel) & test %in% names(tl)]
st[, `:=`(se = abs(gbar / t), test_lab = factor(tl[test], levels = tl),
          fit = factor(sprintf("K = %d", K), levels = sprintf("K = %d", c(K_ref, K_sel))),
          stratum = gsub("_", " ", stratum))]
st[, stratum := factor(stratum, levels = rev(unique(stratum)))]
p4a <- ggplot(st, aes(x = gbar, y = stratum, shape = fit)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbar(aes(xmin = gbar - 2 * se, xmax = gbar + 2 * se), orientation = "y",
                width = 0.3, colour = "grey45", position = position_dodge(width = 0.6)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  facet_wrap(~test_lab, ncol = 1, scales = "free_y") +
  labs(x = expression(bar(g)[b] ~ "(probability mass per word)"), y = NULL,
       title = "Per-stratum mean moments", subtitle = "Whiskers: +/- 2 SE") +
  theme_paper()
tk <- rs$tests$tests[test %in% names(tl)]
tk[, `:=`(sig = pval < 0.05, test_lab = factor(tl[test], levels = tl))]
p4b <- ggplot(tk, aes(x = factor(K), y = gbar_absmax, group = test_lab, linetype = test_lab)) +
  geom_line(linewidth = 0.6) + geom_point(aes(shape = sig), size = 2.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "p < 0.05",
                     labels = c(`TRUE` = "reject", `FALSE` = "not rejected")) +
  scale_y_log10() +
  guides(linetype = guide_legend(nrow = 3), shape = guide_legend(nrow = 2)) +
  labs(x = "Number of topics (K)", y = expression(max[b] ~ "|" * bar(g)[b] * "|" ~ "(log scale)"),
       linetype = NULL, title = "Largest moment at the tested fits",
       subtitle = "Filled points reject at the 5% level") +
  theme_paper()
save_fig(p4a + p4b + plot_layout(widths = c(1, 1.15)), "R4_moment_diagnostics", "MDNA",
         width = 11, height = 5.6, label = L_MDNA)

# ------------------------------- E1: F2 ---------------------------------------------
e1 <- qs_read(p_data("E1", "e1_results_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"))
gm_labs <- function(dt) dt[, `:=`(
  eval_lab = droplevels(factor(eval_label(eval), levels = EVAL_LEVELS)),
  metric_lab = factor(metric_label(metric), levels = METRIC_LEVELS))]
g1 <- gm_labs(copy(e1$gains)[replicate == 1L])
g1[, M := .N, by = c("metric", "eval")]
g1[, ub_simul := delta_mean + qnorm(1 - ALPHA / M) * se]
ex <- rbindlist(lapply(c("ho_reconstruction", "ho_completion"), function(mn)
  select_k_all_rules(e1$doc_rep1[eval == mn], 0.01, ALPHA)[, eval := mn]), use.names = TRUE)
ex <- gm_labs(ex[rule != "adjacent_pointwise" & !is.na(K_hat)])
ex[, rule_lab := factor(ifelse(rule == "total_gain", "total gain", "simultaneous adjacent"),
                        levels = c("simultaneous adjacent", "total gain"))]
p2 <- ggplot(g1, aes(x = K, y = delta_mean)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_hline(yintercept = c(0.01, 0.005), linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lwr, ymax = ub_simul), width = 0, colour = "grey65", linewidth = 1.1) +
  geom_errorbar(aes(ymin = lwr, ymax = ub_onesided), width = 0.3, colour = "grey25") +
  geom_point(size = 1.5) +
  geom_vline(data = ex, aes(xintercept = K_hat, linetype = rule_lab), colour = "grey20") +
  scale_linetype_manual(values = c("simultaneous adjacent" = "dotdash", "total gain" = "solid"),
                        name = "Selection (eps = 0.01)") +
  facet_grid(eval_lab ~ metric_lab, scales = "free_y") +
  labs(x = "Number of topics (K); gain to the next grid point", y = expression(Delta * R^2),
       title = "Held-out adjacent gains with upper bounds (one training seed)",
       subtitle = "Thin whisker: pointwise one-sided bound; thick grey: simultaneous bound; dashed: eps = 0.01, 0.005") +
  theme_paper()
save_fig(p2, "F2_adjacent_gains", "E1", width = 11, height = 6, label = L_SIM)

# ------------------------------- E4: F4, F5 -----------------------------------------
e4 <- qs_read(p_data("E4", "e4_results_E4_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.qs2"))
alt_lab <- c(contamination = "Contamination (refit)",
             contamination_eval = "Contamination (evaluation only)",
             burstiness = "Burstiness (refit); tau", drift = "Vocabulary drift (evaluation only)",
             ctm = "Exchangeable logistic-normal (refit); rho",
             group_vocab = "Group vocabulary (refit)")
T12 <- c(T1_freq_contrast = "Test 1", T2_freq_strata = "Test 2")
K_true <- e4$config$K_true
pw <- e4$power_tab[K == K_true & test %in% names(T12)]
pw[, `:=`(alt_lab = factor(alt_lab[alt], levels = alt_lab), test_lab = T12[test])]
sz <- e4$size_tab[K == K_true & test %in% names(T12)][, test_lab := T12[test]]
p4 <- ggplot(pw, aes(x = strength, y = power, colour = test_lab, group = test_lab)) +
  geom_hline(data = sz, aes(yintercept = size_raw, colour = test_lab), linetype = "dotted") +
  geom_hline(yintercept = 0.05, colour = "grey70") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  facet_wrap(~alt_lab, scales = "free_x") +
  scale_colour_grey(start = 0, end = 0.55, name = NULL) +
  labs(x = "Strength of the alternative", y = "Raw rejection rate (5% level)",
       title = sprintf("Rejection of the moment tests at K = %d", K_true),
       subtitle = "Dotted: baseline rejection under LDA data (not test size); grey: nominal 5%. Test 3 omitted pending recomputation.") +
  theme_paper()
save_fig(p4, "F4_power_curves", "E4", width = 10.4, height = 5.4, label = L_SIM)

r2a <- e4$r2_tab[K == K_true]
pwm <- pw[, .(power = mean(power)), by = c("alt", "strength")]
both <- merge(r2a, pwm, by = c("alt", "strength"))
b_long <- melt(both, id.vars = c("alt", "strength"), measure.vars = c("r2_macro", "power"),
               variable.name = "what")
b_long[, `:=`(what = fifelse(what == "r2_macro", "Held-out Macro R2 (Deviance)",
                             "Raw rejection rate, mean of Tests 1-2"),
              alt_lab = factor(alt_lab[alt], levels = alt_lab))]
p5 <- ggplot(b_long, aes(x = strength, y = value, linetype = what)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  facet_wrap(~alt_lab, scales = "free_x") +
  labs(x = "Strength of the alternative", y = NULL, linetype = NULL,
       title = "Overall fit and moment-test rejection under the alternatives") +
  theme_paper()
save_fig(p5, "F5_fit_vs_tests", "E4", width = 10.4, height = 5.4, label = L_SIM)
log_msg("=== revision figures written under Results/Figures/{%s,%s} ===", L_MDNA, L_SIM)
