# =============================================================================
# postprocess_e2_gap.R -- diagnose the Micro-Macro gap CI undercoverage.
# Pure post-processing of a saved E2 results object; no simulation, no fits.
#
# Table S3 shows systematic 0.90-0.91 coverage of the delta-method gap
# interval at J_ev = 250 (ten of ten K cells below nominal). This script
# separates the two candidate causes:
#   (a) plug-in SE optimism: compare the mean delta-method SE against the
#       empirical SD of the gap estimates WITHIN (K, J_ev, train_seed) --
#       the correct comparison, since the estimand gap_true is conditional
#       on the training fit;
#   (b) bias: mean(gap - gap_true) relative to the empirical SD.
#
# Usage: Rscript Code/postprocess_e2_gap.R [file=Data/E2/e2_results_<TAG>.qs2]
# Output: Results/csv/e2_gap_diagnosis_<TAG>.csv + console conclusion.
# =============================================================================

suppressMessages({library(qs2); library(data.table)})
source(here::here("Code", "R", "source_all.R"))

.args <- commandArgs(trailingOnly = TRUE)
f <- sub("^file=", "", grep("^file=", .args, value = TRUE))
if (length(f) != 1L || !file.exists(f)) stop("Provide file=<explicit existing result.qs2>")
# outputs inherit the input object's revision suffix (e.g. "_rev1"), so
# post-processing a re-scored object never overwrites the original csv files
SFX <- if (grepl("_rev[0-9A-Za-z]*\\.qs2$", f))
  sub("^.*(_rev[0-9A-Za-z]*)\\.qs2$", "\\1", f) else ""
cat(sprintf("post-processing: %s\n", basename(f)))
x <- qs_read(f)
cfg <- x$config
g <- x$gap
stopifnot(all(c("K", "gap", "se", "cover", "gap_true", "rep", "train_seed",
                "J_ev") %in% names(g)))

within <- g[, .(var_within = var(gap), mean_se = mean(se),
                bias = mean(gap - gap_true), n_rep = .N),
            by = .(K, J_ev, train_seed)]
diag <- within[, .(sd_emp    = sqrt(mean(var_within)),
                   mean_se   = mean(mean_se),
                   bias_mean = mean(bias)),
               by = .(K, J_ev)]
diag[, se_ratio := mean_se / sd_emp]            # < 1 => SE optimistic
diag[, bias_over_sd := bias_mean / sd_emp]

zt <- g[, .(cover_emp = mean(cover),
            cover_z   = mean(abs((gap - gap_true) / se) <= qnorm(0.975)),
            share_z_hi = mean((gap - gap_true) / se >  qnorm(0.975)),
            share_z_lo = mean((gap - gap_true) / se < -qnorm(0.975))),
        by = .(K, J_ev)]
diag <- merge(diag, zt, by = c("K", "J_ev"))[order(J_ev, K)]

write_result(diag, paste0("e2_gap_diagnosis", SFX), cfg)
cat("\n== gap CI diagnosis by (J_ev, K) ==\n")
print(diag, nrows = 60L, digits = 3)

s <- diag[, .(cover = mean(cover_emp), se_ratio = mean(se_ratio),
              bias_over_sd = mean(bias_over_sd),
              hi = mean(share_z_hi), lo = mean(share_z_lo)), by = J_ev]
cat("\n== summary by J_ev ==\n"); print(s, digits = 3)
cat("\nREADING GUIDE:\n")
cat(" - se_ratio < 1 at J_ev = 250 while ~1 elsewhere  => plug-in SE optimism.\n")
cat(" - |bias_over_sd| large                            => bias, not SE.\n")
cat(" - hi vs lo asymmetry                              => one-sided miss\n")
cat("   (e.g. skewness of the ratio estimator at small J_ev).\n")
log_msg("postprocess_e2_gap done for %s", basename(f))
