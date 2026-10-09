# =============================================================================
# make_revision_tex.R -- LaTeX table BODIES for the Sept-2026 revision, written
# from the gated csv outputs of postprocess_revision.R / rescore_revision.R so
# that no number in the revised supplement is transcribed by hand.
# Output: Results/tex/rev1_<table>.tex (rows only; captions live in the paper).
# Usage:  Rscript Code/make_revision_tex.R [suffix=_rev1]
# =============================================================================
suppressMessages({library(data.table); library(qs2)})
source(here::here("Code", "R", "source_all.R"))
.args <- commandArgs(trailingOnly = TRUE)
SFX <- sub("^suffix=", "", grep("^suffix=", .args, value = TRUE))
if (!length(SFX)) SFX <- "_rev1"
rd <- function(name) fread(p_results("csv", sprintf("%s%s.csv", name, SFX)))
put <- function(lines, name) {
  f <- p_results("tex", sprintf("rev1_%s.tex", name))
  writeLines(lines, f); log_msg("tex %s (%d rows)", basename(f), length(lines))
}
f4 <- function(x) ifelse(x < 0, sprintf("$-%.4f$", abs(x)), sprintf("$%.4f$", x))
sci <- function(x) { e <- floor(log10(abs(x))); m <- x / 10^e
  ifelse(abs(x) >= 1e-3, sprintf("$%.*f$", ifelse(abs(x) >= 0.1, 2, 3), x),
         sprintf("$%.1f\\times10^{%d}$", m, e)) }
kh <- function(k) ifelse(is.na(k), "none", as.character(k))

# --- S13: gain profile with the three bounds ------------------------------------
p <- rd("mdna_gain_profile")
w <- dcast(p, K + K_to ~ protocol, value.var = c("delta_mean", "ub_pointwise",
                                                 "ub_adjacent_simul", "total_gain_ub"))
put(w[order(K), sprintf("$%d\\to%d$ & %s & %s & %s & %s & %s & %s & %s & %s \\\\", K, K_to,
    f4(delta_mean_com), f4(ub_pointwise_com), f4(ub_adjacent_simul_com), f4(total_gain_ub_com),
    f4(delta_mean_rec), f4(ub_pointwise_rec), f4(ub_adjacent_simul_rec), f4(total_gain_ub_rec))],
    "S13_gain_profile")

# --- selections by family (completion + reconstruction) -------------------------
s <- rd("mdna_selection_rules")
lab_rule <- c(total_gain = "Total gain, simultaneous",
              adjacent_simultaneous = "Adjacent gain, simultaneous",
              adjacent_pointwise = "Adjacent gain, pointwise")
lab_fam <- c(dev = "Deviance", chisq = "Pearson", se = "Squared-Error")
sw <- dcast(s, metric + rule ~ protocol + eps, value.var = "K_hat")
sw[, `:=`(o1 = match(metric, names(lab_fam)), o2 = match(rule, names(lab_rule)))]
put(sw[order(o1, o2), sprintf("%s & %s & %s & %s & %s & %s \\\\", lab_fam[metric], lab_rule[rule],
    kh(com_0.01), kh(com_0.005), kh(rec_0.01), kh(rec_0.005))], "S13b_selection_by_family")

# --- S16: stratum-level moments, reference block and selection block -------------
st <- rd("mdna_moment_strata")[test %in% c("T2_freq_strata", "T3_fit_strata")]
st[, lab := sub("_vs_", "$ vs.\\\\ $", stratum)]
st[, lab := sprintf("$%s$", gsub("([fs])([0-9])", "\\1_\\2", lab))]
blk <- function(Ks, name) {
  w <- dcast(st[K %in% Ks], test + stratum + lab ~ K, value.var = c("gbar", "t"))
  w <- w[order(test, stratum)]
  cols <- unlist(lapply(Ks, function(k) c(sprintf("gbar_%d", k), sprintf("t_%d", k))))
  rows <- apply(w, 1, function(r) {
    v <- vapply(seq_along(cols), function(i) {
      x <- as.numeric(r[[cols[i]]])
      if (i %% 2 == 1) sprintf("$%s%.1f$", ifelse(x < 0, "-", ""), abs(x * 1e7))
      else sprintf("$%s%.1f$", ifelse(x < 0, "-", ""), abs(x)) }, "")
    sprintf("%s & %s & %s \\\\", ifelse(r[["test"]] == "T2_freq_strata", "Test 2", "Test 3"),
            r[["lab"]], paste(v, collapse = " & "))
  })
  put(rows, name)
}
blk(c(40, 50, 60), "S16_strata_reference"); blk(c(170, 180, 190), "S16_strata_selection")

# --- joint tests, both covariances, all six K -------------------------------------
tt <- merge(rd("mdna_moment_tests")[test != "T3_bug_minK_strata"],
            rd("mdna_moment_tests_cluster")[test != "T3_bug_minK_strata",
                                            .(K, test, stat_cluster, pval_cluster)],
            by = c("K", "test"))
tl <- c(T1_freq_contrast = "Test 1", T2_freq_strata = "Test 2", T3_fit_strata = "Test 3")
put(tt[order(test, K), sprintf("%s & %d & %.1f & %s & %.1f & %s \\\\", tl[test], K, stat,
    sci(pval), stat_cluster, sci(pval_cluster))], "S16b_joint_tests")

# --- residual mass by group -----------------------------------------------------
m <- rd("mdna_residual_mass")[K %in% c(50, 180)]
mw <- dcast(m, partition + stratum_lab + words ~ K, value.var = "mass_pp")
put(mw[order(partition, stratum_lab), sprintf("%s & %s & %d & $%+.4f$ & $%+.4f$ \\\\",
    sub(" \\(.*", "", partition), stratum_lab, words, `50`, `180`)], "S16c_residual_mass")

# --- S17: residual vocabulary (unchanged) + worst fit under the Poisson null -------
x <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016.qs2"))
up <- head(x$resid_words[order(-resid_mean)], 8L)
wf <- head(rd("mdna_wordlevel_worst_fit_paper_null"), 8L)
put(sprintf("%s & %.1f & %.2f & %s (%d) & $%.2f$ \\\\", up$word, up$resid_mean * 1e4,
            100 * up$share_train, wf$word, wf$doc_freq, wf$r2_paper_null), "S17_vocabulary")

# --- word-level curve: Micro under the Poisson null for every K ---------------------
mk <- rd("mdna_wordlevel_micro_over_K")
put(mk[K %in% c(10, 50, 100, 150, 200), sprintf("%d & %.2f & %.2f \\\\", K,
    w_micro_paper_null, w_micro_package)], "S15_word_micro")

# --- delta panel, design panel, firm panel (Table S19) -----------------------------
d <- rd("mdna_delta_sensitivity")[protocol == "com"][order(delta)]
put(d[, sprintf("%s & %d & %.3f & %.3f & %s & %s & %s & %s & %s & %s \\\\",
    ifelse(delta == 0, "0 (strictly positive)", format(delta)), retained, macro_K50, micro_K50,
    kh(total_gain_0.01), kh(total_gain_0.005), kh(adjacent_simultaneous_0.01),
    kh(adjacent_simultaneous_0.005), kh(adjacent_pointwise_0.01), kh(adjacent_pointwise_0.005))],
    "S19_delta_panel")
g <- rd("mdna_selection_design_sensitivity")
gw <- dcast(g, design + n_models_in_support + max_abs_dR2 + rule ~ protocol + eps,
            value.var = "K_hat")
gw[, o := match(rule, names(lab_rule))]
put(gw[order(design, o), sprintf("%s & %s & %s & %s & %s & %s \\\\", design, lab_rule[rule],
    kh(com_0.01), kh(com_0.005), kh(rec_0.01), kh(rec_0.005))], "S19_design_panel")
fs <- rd("mdna_selection_firm_sensitivity")
fw <- dcast(fs, sample + rule ~ protocol + eps, value.var = "K_hat")
fw[, o := match(rule, names(lab_rule))]
put(fw[order(sample, o), sprintf("%s & %s & %s & %s & %s & %s \\\\", sample, lab_rule[rule],
    kh(com_0.01), kh(com_0.005), kh(rec_0.01), kh(rec_0.005))], "S19_firm_panel")
log_msg("=== make_revision_tex complete ===")
