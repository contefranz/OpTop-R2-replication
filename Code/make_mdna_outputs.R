# =============================================================================
# make_mdna_outputs.R -- the MINIMAL paper exhibit set for the MD&A study.
# Reads Data/MDNA/mdna_results_MDNA_<y1>_<y2>[_n<sample>].qs2 (run_mdna.R).
#
# Paper set (3 figures + 3 tables):
#   R1  figure: fit curves -- in-sample / reconstruction / completion,
#       three family facets, K-hat marked (the consistency picture)
#   R1b figure: Micro - Macro gap over K (delta-method CI, reconstruction)
#   R1c figure: adjacent gains over K x held-out target, with one-sided upper
#       bounds and the eps-rule reference lines
#   R4  figure: moment diagnostics -- per-stratum mean moments at K-hat (where
#       mass is misallocated) + effect sizes across the tested grid (does more
#       K repair the violations; needs run_mdna tests_K=all for the full grid)
#   A3  figure: doc-level fit vs length at K-hat (appendix)
# families=dev (default) renders Deviance-only figures (family agreement is
# certified by the R1 selection table); families=all restores 3-family facets.
#   R1  table : selection -- eps-rules (families x targets x eps) + comparators
#   R2  table : consistency battery at K-hat -- Micro/Macro + Prop-2 CI per
#       family & protocol; Micro-Macro gap + Prop-1(iii) channels; null-floor
#       (boilerplate) share; Section-4 moment tests (stat, p, gbar)
#   R3  table : vocabulary -- top over-observed words (signed residual) and
#       worst-fit words (word-level R2)
# Appendix extras (generated, clearly named):
#   A1  figure: word-level dual (w-Micro/w-Macro over K)
#   A2  table : boilerplate share by FF12 industry
#
# Usage: Rscript Code/make_mdna_outputs.R [y1=2015] [y2=2016] [sample_n=0]
# =============================================================================

suppressMessages({
  library(data.table); library(ggplot2); library(patchwork); library(qs2)
  library(tinytable); library(writexl)
})
source(here::here("Code", "R", "source_all.R"))

.args <- commandArgs(trailingOnly = TRUE)
P <- list(y1 = 2015L, y2 = 2016L, sample_n = 0L,
          families = "dev")   # "dev" = Deviance-only figures; "all" = 3 families
for (a in .args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1L]]
  if (length(kv) != 2L || !kv[1L] %in% names(P))
    stop("unknown argument '", a, "'. Valid: ",
         paste(names(P), collapse = ", "), call. = FALSE)
  P[[kv[1L]]] <- if (kv[1L] == "families") kv[2L] else as.integer(kv[2L])
}
stopifnot(P$families %in% c("dev", "all"))
DEV_ONLY <- P$families == "dev"
fam_note <- if (DEV_ONLY)
  "Deviance Macro Index" else
  "all three discrepancy families"
suffix <- if (P$sample_n > 0L) sprintf("_n%d", P$sample_n) else ""
TAG <- sprintf("MDNA_%d_%d%s", P$y1, P$y2, suffix)
RUN_LABEL <- sprintf("mdna_%d_%d%s", P$y1, P$y2, suffix)
f <- p_data("MDNA", sprintf("mdna_results_%s.qs2", TAG))
stopifnot(file.exists(f))
x <- cache_get(f)
need <- c("summary", "khat", "sel_comparators", "battery_ci", "gap", "decomp",
          "bp_doc", "tests_ho", "resid_words", "wstar", "wcurve", "boilerplate")
miss <- need[vapply(need, function(nm) is.null(x[[nm]]), logical(1L))]
if (length(miss))
  stop("results file ", basename(f), " predates the current pooled-design ",
       "run_mdna.R (missing fields: ", paste(miss, collapse = ", "),
       "). Re-run prep_mdna.R and run_mdna.R with the same y1/y2/sample_n.",
       call. = FALSE)
K_hat <- x$K_hat
log_msg("MDNA outputs for %s (pooled fiscal %d+%d, K-hat = %d)",
        TAG, P$y1, P$y2, K_hat)

.tabs <- list()
emit <- function(dt, name) {
  .tabs[[name]] <<- as.data.frame(dt)
  fwrite(dt, p_results("csv", sprintf("%s_%s.csv", name, RUN_LABEL)))
  tryCatch(
    save_tt(tt(as.data.frame(dt)),
            p_results("tex", sprintf("%s_%s.tex", name, RUN_LABEL)),
            overwrite = TRUE),
    error = function(e) log_msg("tex skipped for %s: %s", name,
                                conditionMessage(e)))
  log_msg("table %s written", name)
}
eval_lab_m <- c(ins = "In-sample", rec = "Held-out (reconstruction)",
                com = "Held-out (completion)")

# =========================== R1: fit curves (figure) ===========================
s <- copy(x$summary)
if (DEV_ONLY) s <- s[metric == "dev"]
s[, eval_lab := factor(eval_lab_m[eval], levels = EVAL_LEVELS)]
s[, metric_lab := factor(metric_label(metric), levels = METRIC_LEVELS)]
p1 <- ggplot(s, aes(x = K, y = r2_macro, linetype = eval_lab)) +
  geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.5) +
  geom_vline(xintercept = K_hat, linetype = "dotdash", colour = "grey40") +
  {if (!DEV_ONLY) facet_wrap(~metric_lab)} +  # fixed scales across families
  scale_eval_linetype() + scale_eval_shape() +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Number of topics (K)", y = expression(R["Macro"]^2),
       title = sprintf("MD&A corpus (fiscal %d-%d): fit over K", P$y1, P$y2),
       subtitle = sprintf(
         "In-sample and held-out protocols; K-hat = %d (completion, eps = 0.01)",
         K_hat)) +
  theme_paper()
save_fig(p1, "R1_fit_curves", "MDNA",
         width = if (DEV_ONLY) 7.4 else 10.2, height = 4.4)

# ================= R1b: Micro - Macro gap over K (figure) ======================
g <- copy(x$gap)
if (DEV_ONLY) g <- g[metric == "dev"]
g[, eval_lab := factor(eval_lab_m[eval], levels = EVAL_LEVELS)]
g[, metric_lab := factor(metric_label(metric), levels = METRIC_LEVELS)]
p1b <- ggplot(g, aes(x = K, y = gap, linetype = eval_lab)) +
  geom_ribbon(data = g[eval == "rec"],
              aes(x = K, ymin = lwr, ymax = upr),
              inherit.aes = FALSE,     # keep the ribbon out of the linetype
              fill = "grey85", alpha = 0.6, colour = NA) +  # guide so the
  # shape and linetype guides stay merged into one "Evaluation" legend
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) + geom_point(aes(shape = eval_lab), size = 1.5) +
  geom_vline(xintercept = K_hat, linetype = "dotdash", colour = "grey40") +
  {if (!DEV_ONLY) facet_wrap(~metric_lab)} +
  scale_eval_linetype() + scale_eval_shape() +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Number of topics (K)", y = "Micro - Macro gap",
       title = "Aggregation heterogeneity: Micro - Macro gap over K",
       subtitle = sprintf(
         "%s; ribbon: delta-method 95%% CI on the reconstruction target",
         fam_note)) +
  theme_paper()
save_fig(p1b, "R1b_micro_macro_gap", "MDNA",
         width = if (DEV_ONLY) 7.4 else 10.2, height = 4.4)

# ============ R1c: adjacent gains + eps rule (figure) ==========================
ga <- copy(x$gains)
if (DEV_ONLY) ga <- ga[metric == "dev"]
ga[, eval_lab := factor(eval_lab_m[eval], levels = EVAL_LEVELS)]
ga[, metric_lab := factor(metric_label(metric), levels = METRIC_LEVELS)]
p1c <- ggplot(ga, aes(x = K_next, y = delta_mean)) +
  geom_hline(yintercept = c(0.01, 0.005), linetype = "dashed",
             colour = "grey55") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_linerange(aes(ymin = delta_mean, ymax = ub_onesided), colour = "grey45") +
  geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
  geom_vline(xintercept = K_hat, linetype = "dotdash", colour = "grey40") +
  {if (DEV_ONLY) facet_wrap(~eval_lab, ncol = 2)
   else facet_grid(eval_lab ~ metric_lab)} +
  labs(x = "Number of topics (K)",
       y = expression(Delta ~ R["Macro"]^2 ~ "(gain from previous grid point)"),
       title = "Adjacent held-out gains and the eps-adequacy rule",
       subtitle = sprintf(
         "%s; whiskers: one-sided upper bounds; dashed lines: eps = 0.01 and 0.005",
         fam_note)) +
  theme_paper()
save_fig(p1c, "R1c_adjacent_gains", "MDNA",
         width = if (DEV_ONLY) 9.4 else 10.2,
         height = if (DEV_ONLY) 4.2 else 6.2)

# =========================== R1: selection (table) =============================
sel <- rbindlist(list(
  x$khat[, .(rule = sprintf("eps-rule %s (%s, %s)",
                            sub("eps_", "", rule), metric, eval), K_hat)],
  x$sel_comparators[, .(rule, K_hat)]), fill = TRUE)
emit(sel, "R1_selection")

# ================== R2: consistency battery at K-hat (table) ===================
bat1 <- x$battery_ci[, .(block = "fit", eval = eval_lab_m[eval], metric,
                         value = round(r2_macro, 4),
                         ci = sprintf("[%.4f, %.4f]", lwr, upr))]
mic <- x$summary[K == K_hat, .(block = "fit_micro", eval = eval_lab_m[eval],
                               metric, value = round(r2_micro, 4), ci = "")]
gapK <- x$gap[K == K_hat, .(block = "gap", eval = eval_lab_m[eval], metric,
                            value = round(gap, 4),
                            ci = sprintf("[%.4f, %.4f]", lwr, upr))]
dec <- x$decomp[K == K_hat,
                .(block = "gap_channels", eval = "Held-out (reconstruction)",
                  metric = c("length", "atypicality", "interaction"),
                  value = round(c(ch_length, ch_atypicality, ch_interaction), 4),
                  ci = "")]
bp <- data.table(block = "boilerplate", eval = "Held-out (completion)",
                 metric = "excluded_share",
                 value = round(mean(x$bp_doc$excluded), 4), ci = "")
tst <- x$tests_ho[K == K_hat,
                  .(block = "moment_tests", eval = "Held-out", metric = test,
                    value = round(stat, 1),
                    ci = sprintf("p = %.2g; gbar = %.2g", pval, gbar_absmax))]
emit(rbindlist(list(bat1, mic, gapK, dec, bp, tst), fill = TRUE), "R2_battery")

# ======================= R3: vocabulary (table) =================================
over <- head(x$resid_words, 15L)[, .(criterion = "over-observed (residual)",
                                     word, value = signif(resid_mean, 3))]
worst <- head(x$wstar[keep == TRUE & !is.na(r2_word)][order(r2_word)],
              15L)[, .(criterion = "worst fit (word R2)",
                                word = word_id, value = signif(r2_word, 3))]
emit(rbindlist(list(over, worst)), "R3_vocabulary")

# ============== R4: moment diagnostics (figure, two panels) =====================
# (a) WHERE the model misallocates mass: per-stratum mean moments at K-hat,
#     whiskers +/- 2 SE (SE recovered from the reported t-statistic).
st <- x$strata_ho[K == K_hat]
st[, se := fifelse(abs(t) > 1e-8, abs(gbar / t), NA_real_)]
st[, stratum := factor(stratum, levels = rev(unique(stratum)))]
p4a <- ggplot(st, aes(x = gbar, y = stratum)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_errorbar(aes(xmin = gbar - 2 * se, xmax = gbar + 2 * se),
                orientation = "y", width = 0.25, colour = "grey45",
                na.rm = TRUE) +
  geom_point(size = 1.8) +
  facet_wrap(~test, ncol = 1, scales = "free_y") +
  labs(x = expression(bar(g)[b] ~ "(probability-mass units)"), y = NULL,
       title = sprintf("Per-stratum mean moments at K = %d", K_hat),
       subtitle = "Whiskers: +/- 2 SE") +
  theme_paper()
# (b) does more K repair the violations? effect size over the tested grid,
#     filled points = rejection at 5%.
tk <- copy(x$tests_ho)
tk[, sig := pval < 0.05]
p4b <- ggplot(tk, aes(x = K, y = gbar_absmax, linetype = test)) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 2.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     name = "p < 0.05", labels = c(`TRUE` = "reject",
                                                   `FALSE` = "not rejected")) +
  scale_y_log10() +
  geom_vline(xintercept = K_hat, linetype = "dotdash", colour = "grey40") +
  labs(x = "Number of topics (K)",
       y = expression(max[b] ~ "|" * bar(g)[b] * "|" ~ "(log scale)"),
       linetype = "Test",
       title = "Moment violations across the grid",
       subtitle = "Effect sizes; filled points reject at the 5% level") +
  theme_paper()
save_fig(p4a + p4b + plot_layout(widths = c(1, 1.3)),
         "R4_moment_diagnostics", "MDNA", width = 11, height = 5.2)

# ===================== appendix extras (clearly named) ==========================
wc <- melt(x$wcurve[, .(K, w_Micro = r2_micro_word, w_Macro = r2_macro_word)],
           id.vars = "K", variable.name = "agg", value.name = "r2")
pa1 <- ggplot(wc, aes(x = K, y = r2, linetype = agg)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  geom_vline(xintercept = K_hat, linetype = "dotdash", colour = "grey40") +
  scale_linetype_manual(values = c(w_Micro = "solid", w_Macro = "dotted"),
                        name = "Aggregation") +
  labs(x = "Number of topics (K)", y = expression(R["Dev,w"]^2),
       title = "Word-level dual perspective (appendix)") +
  theme_paper()
save_fig(pa1, "A1_word_micro_macro", "MDNA", width = 7.4)
emit(x$boilerplate, "A2_boilerplate_by_industry")

# A3: the mechanism behind R1b -- doc-level fit vs length at K-hat
dl <- x$doc_rec[metric == "dev" & K == K_hat & !is.na(r2_doc)]
pa3 <- ggplot(dl, aes(x = L, y = r2_doc)) +
  geom_point(size = 0.8, alpha = 0.35, colour = "grey55") +
  stat_summary_bin(fun = median, bins = 12, geom = "line",
                   colour = "black", linewidth = 0.8) +
  scale_x_log10() +
  labs(x = "Document length (tokens, log scale)",
       y = bquote(R["Dev,j"]^2 ~ "at K =" ~ .(K_hat)),
       title = "Document-level fit vs length (appendix)",
       subtitle = "Line: binned medians -- the length channel behind the Micro-Macro gap") +
  theme_paper()
save_fig(pa3, "A3_fit_vs_length", "MDNA", width = 7.4)

writexl::write_xlsx(.tabs, p_results("xlsx", sprintf("MDNA_tables_%d_%d%s.xlsx",
                                                     P$y1, P$y2, suffix)))
log_msg("workbook: MDNA_tables_%d_%d%s.xlsx (%d sheets)", P$y1, P$y2, suffix,
        length(.tabs))
log_msg("=== MDNA outputs complete (label %s) ===", RUN_LABEL)
