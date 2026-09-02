# =============================================================================
# make_tables.R -- all Section 5 tables (T1-T8) in three formats:
#   * LaTeX fragments (tinytable)      -> Results/tex/<name>_<profile>.tex
#   * CSV copies                       -> Results/csv/<name>_<profile>.csv
#   * ONE Excel workbook per profile   -> Results/xlsx/Section5_tables_<profile>.xlsx
#     (one sheet per table, unrounded values)
# Usage:  Rscript Code/make_tables.R [smoke|pilot|full]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
source(here::here("Code", "R", "source_all.R"))
source(here::here("Code", "config", "configs.R"))
cli <- parse_cli(args)
PROFILE <- cli$profile
sel_exp <- setdiff(toupper(strsplit(cli$overrides$exp %||% "", ",")[[1]]), "")
cli$overrides$exp <- NULL
RUN_LABEL <- run_label(list(profile = PROFILE), cli$overrides, cli$label)
log_msg("tables for run label '%s'%s", RUN_LABEL,
        if (length(sel_exp)) sprintf(" (experiments: %s)",
                                     paste(sel_exp, collapse = ",")) else "")
library(tinytable)

.tab_registry <- list()

load_results <- function(exp) {
  if (length(sel_exp) && !toupper(exp) %in% sel_exp) return(NULL)
  cfg <- apply_overrides(get_config(exp, PROFILE), cli$overrides, strict = FALSE)
  f <- p_data(exp, sprintf("%s_results_%s.qs2", tolower(exp), run_tag(cfg)))
  if (!file.exists(f)) { log_msg("skip %s tables", exp); return(NULL) }
  cache_get(f)
}

emit <- function(dt, name, caption = NULL, digits = 4) {
  .tab_registry[[name]] <<- as.data.frame(dt)
  fwrite(dt, p_results("csv", sprintf("%s_%s.csv", name, RUN_LABEL)))
  num <- names(dt)[vapply(dt, is.numeric, TRUE)]
  dt_r <- copy(dt)[, (num) := lapply(.SD, function(x)
    ifelse(abs(x) < 5e-5 & x != 0 & abs(x) > 0, signif(x, 2), round(x, digits))),
    .SDcols = num]
  tt(as.data.frame(dt_r), caption = caption) |>
    save_tt(p_results("tex", sprintf("%s_%s.tex", name, RUN_LABEL)),
            overwrite = TRUE)
  log_msg("table %s written", name)
}

e1 <- load_results("E1"); e2 <- load_results("E2"); e3 <- load_results("E3")
e4 <- load_results("E4"); e5 <- load_results("E5"); e6 <- load_results("E6")

# --- T1: K-hat selection distribution (E1) -----------------------------------------
if (!is.null(e1)) {
  K_true <- e1$config$K_true
  kh <- e1$khat[metric %in% c("dev", "chisq", "se", "comparator")]
  kh[, rule_full := fifelse(eval == "comparator", rule,
                            sprintf("%s_%s_%s", rule, metric, eval))]
  t1 <- kh[, .(
    n = .N,
    mode_K = { tb <- table(K_hat)
               if (length(tb)) as.integer(names(tb)[which.max(tb)]) else NA_integer_ },
    pct_correct = mean(K_hat == K_true, na.rm = TRUE),
    pct_within1 = mean(abs(K_hat - K_true) <= 1, na.rm = TRUE),
    mean_K = mean(K_hat, na.rm = TRUE),
    pct_na = mean(is.na(K_hat))
  ), by = rule_full][order(rule_full)]
  emit(t1, "T1_khat_distribution",
       sprintf("Selection distribution across %d replicates (K* = %d)",
               e1$config$S, K_true))

  # T1b: per-K held-out Macro by family and target + the completion share of
  # documents excluded by the 0.14.1 null-discrepancy floor (D_null < c).
  # Direct refresh source for the draft's tab:e1_r2.
  hb <- e1$summary[eval %in% c("ho_reconstruction", "ho_completion"),
                   .(r2_macro = mean(r2_macro)), by = .(K, metric, eval)]
  t1b <- dcast(hb, K ~ eval + metric, value.var = "r2_macro")
  if ("null_excl_share" %in% names(e1$summary)) {
    ex <- e1$summary[eval == "ho_completion" & metric == "dev",
                     .(completion_excl_share = mean(null_excl_share)), by = K]
    t1b <- merge(t1b, ex, by = "K", all.x = TRUE)
  }
  emit(t1b, "T1b_heldout_by_K",
       sprintf(paste("Held-out Macro indices by K and target; completion share",
                     "of documents excluded by the null floor (K* = %d)"),
               K_true))
}

# --- T2: coverage + paired-gain size/power (E2) --------------------------------------
if (!is.null(e2)) {
  emit(dcast(e2$cover_tab, K ~ J_ev, value.var = "coverage"),
       "T2a_coverage", "Empirical coverage of the 95% held-out Macro CI")
  emit(dcast(e2$size_tab, K ~ J_ev, value.var = "size_centered"),
       "T2b_gain_test_size",
       "Size of the paired adjacent-gain test (centered at conditional truth)")
  emit(dcast(e2$power_tab, K ~ J_ev, value.var = "power_pos"),
       "T2c_gain_test_power", "Power against zero gain (K < K*)")
  emit(dcast(e2$gap_tab, K ~ J_ev, value.var = "coverage"),
       "T2d_gap_coverage", "Coverage of the delta-method gap CI (Remark 6)")
  emit(dcast(e2$khat_tab, eps + K_hat ~ J_ev, value.var = "N", fill = 0L),
       "T2e_khat_by_Jev", "Selection distribution of K-hat by J_ev")
}

# --- T3/T4/T5: moment tests (E4) ------------------------------------------------------
if (!is.null(e4)) {
  emit(e4$size_tab, "T3_moment_size",
       "Tests 1-3 under the correct DGP: raw (H0: mu = 0) and centered-at-conditional-truth rejection rates (nominal 5%), with mean effect sizes")
  if (nrow(e4$power)) {
    emit(dcast(e4$power_tab, alt + strength ~ test, value.var = "power",
               subset = .(K == e4$config$K_true)),
         "T4_moment_power",
         sprintf("Power at K = %d by alternative and strength", e4$config$K_true))
    emit(e4$r2_tab[K == e4$config$K_true],
         "T4b_r2_under_misspec", "Held-out fit under misspecification")
  }
  if (!is.null(e4$word_resid)) {
    wres <- e4$word_resid[order(-resid_mean)]
    n_pl <- sum(wres$planted)
    t5 <- data.table(
      criterion = "mean held-out residual (desc)",
      n_planted = n_pl,
      precision_at_n = mean(head(wres, n_pl)$planted),
      recall_top50 = sum(head(wres, 50)$planted) / n_pl,
      median_rank_planted = median(which(wres$planted)),
      W = nrow(wres)
    )
    if (!is.null(e4$word_rank)) {
      wr <- e4$word_rank[keep == TRUE][order(r2_word)]
      t5 <- rbind(t5, data.table(
        criterion = "word-level R2 (asc, filtered)",
        n_planted = sum(wr$planted),
        precision_at_n = mean(head(wr, sum(wr$planted))$planted),
        recall_top50 = sum(head(wr, 50)$planted) / max(1, sum(wr$planted)),
        median_rank_planted = median(which(wr$planted)),
        W = nrow(wr)
      ))
    }
    emit(t5, "T5_planted_words",
         "Identification of planted contamination words by two word-level criteria")
  }
  if (!is.null(e4$strata_ill)) {
    emit(e4$strata_ill, "T5b_strata_illustration",
         "Per-stratum moments for one rejected case (BH-adjusted)")
  }
}

# --- T6: gap + decomposition + length stats (E3) --------------------------------------
if (!is.null(e3)) {
  K_true <- e3$config$K_true
  # gap at K* by discrepancy family (held-out); Prop.-1(iii) channels are
  # Deviance-only (the 3-channel identity exists only for the deviance family)
  g6 <- e3$gap[K == K_true & eval == "ho_reconstruction",
               .(gap = mean(gap)), by = .(scenario, metric)]
  g6w <- dcast(g6, scenario ~ metric, value.var = "gap")
  setnames(g6w, c("dev", "chisq", "se"),
           c("gap_Dev", "gap_Pearson", "gap_SE"), skip_absent = TRUE)
  d6 <- e3$decomp[K == K_true & eval == "ho_reconstruction",
                  .(dev_len = mean(ch_length),
                    dev_atyp = mean(ch_atypicality),
                    dev_inter = mean(ch_interaction)),
                  by = scenario]
  emit(merge(g6w, d6, by = "scenario"),
       "T6a_gap_decomposition",
       sprintf("Held-out Micro-Macro gap at K* = %d by discrepancy family; Deviance Prop.-1(iii) channels", K_true))
  emit(e3$lenstats[, lapply(.SD, mean), by = scenario,
                   .SDcols = c("mean", "sd", "p25", "median", "p75")],
       "T6b_length_stats", "Document-length distributions by scenario")
}

# --- T7: design sensitivity (E5) -------------------------------------------------------
if (!is.null(e5)) {
  emit(e5$grid_sensitivity, "T7a_grid_sensitivity",
       "Reported index for the same K under nested estimation grids")
  emit(e5$minbin, "T7b_minbin",
       "Min-bin document share and excluded probability mass")
}

# --- T9: word-level dual perspective (E6) ----------------------------------------------
if (!is.null(e6)) {
  K_true <- e6$config$K_true
  t9 <- e6$curve[K == K_true & eval == "ho_reconstruction",
                 .(w_Micro = mean(r2_micro_word),
                   w_Macro = mean(r2_macro_word),
                   gap = mean(r2_micro_word - r2_macro_word),
                   n_words = mean(n_words)), by = scenario]
  l2 <- e6$lemma2[, .(lemma2_resid_max = max(abs(resid))), by = scenario]
  t9 <- merge(t9, l2, by = "scenario")
  if (nrow(e6$word_star)) {
    prec <- e6$word_star[keep == TRUE & planted == 1L, .(n_planted = .N),
                         by = scenario]
    top <- e6$word_star[keep == TRUE][order(scenario, r2_word)]
    hit <- top[, .(precision_at_planted = {
      np <- sum(planted); if (np) mean(head(planted, np)) else NA_real_
    }), by = scenario]
    t9 <- merge(t9, hit, by = "scenario", all.x = TRUE)
  }
  emit(t9, "T9_word_micro_macro",
       sprintf("Word-level w-Micro/w-Macro gap at K* = %d, Lemma-2 residual, planted-word precision",
               K_true))
}

# --- T8: runtime and configuration summary ---------------------------------------------
if (!is.null(e1)) {
  t8 <- e1$diagnostics[, .(fits = .N, core_minutes = sum(elapsed_s) / 60,
                           mean_s_per_fit = mean(elapsed_s)), by = method]
  emit(t8, "T8_runtime", sprintf("LDA fitting cost (profile %s)", PROFILE))
}

# --- Excel workbook: one sheet per table ---------------------------------------
if (length(.tab_registry)) {
  sheets <- .tab_registry
  names(sheets) <- substr(names(sheets), 1L, 31L)   # Excel sheet-name limit
  wb <- p_results("xlsx", sprintf("Section5_tables_%s.xlsx", RUN_LABEL))
  writexl::write_xlsx(sheets, wb)
  log_msg("workbook: %s (%d sheets)", basename(wb), length(sheets))
}

log_msg("tables complete (label %s)", RUN_LABEL)
