# =============================================================================
# postprocess_e4.R -- size-adjusted power + contamination-refit seed profile.
# Pure post-processing of a saved E4 results object; no simulation, no fits.
#
# (i)  Size-adjusted power. The conditional null rejects at raw rates far above
#      nominal (Table S7), so raw power under the refit alternatives is not
#      interpretable as discrimination. This script computes empirical 95%
#      critical values from the correct-DGP null statistics ($size$stat),
#      both pooled per (test, K) and seed-matched per (test, K, train_seed)
#      -- the seed-matched version respects the conditional-on-training-fit
#      null (500 null draws per seed) -- and applies them to the alternative
#      statistics ($power$stat).
# (ii) Contamination-refit profile by seed: per-(strength, train_seed) mean
#      held-out R2 at K = K_true, to explain (or expose) the non-monotone
#      aggregate fit across strengths (0.460 -> 0.375 -> 0.502).
#
# Usage: Rscript Code/postprocess_e4.R [file=Data/E4/e4_results_<TAG>.qs2]
#        (defaults to the newest full-profile E4 results object)
# Outputs: Results/csv/e4_power_sizeadj_<TAG>.csv
#          Results/csv/e4_contamrefit_by_seed_<TAG>.csv
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
K_true <- cfg$K_true

# ------------------------- (i) size-adjusted power ----------------------------
size <- x$size; power <- x$power
stopifnot(all(c("test", "stat", "pval", "K", "train_seed") %in% names(size)))
HAS_POWER <- !is.null(power) && nrow(power) > 0L
if (!HAS_POWER) cat("NOTE: no power arm in this object (size-only run);",
                    "skipping the power adjustment and the seed profile.\n")

crit_pool <- size[, .(crit_pool = quantile(stat, 0.95, type = 7)),
                  by = .(test, K)]
crit_seed <- size[, .(crit_seed = quantile(stat, 0.95, type = 7)),
                  by = .(test, K, train_seed)]

if (HAS_POWER) {
pw <- merge(power, crit_pool, by = c("test", "K"))
pw <- merge(pw, crit_seed, by = c("test", "K", "train_seed"))
adj <- pw[, .(power_raw       = mean(pval < 0.05),
              power_adj_pool  = mean(stat > crit_pool),
              power_adj_seed  = mean(stat > crit_seed),
              n = .N),
          by = .(alt, strength, test, K)][order(alt, strength, test, K)]
} # HAS_POWER

# sanity: reproduce the published raw and centred sizes
size_check <- size[, .(size_raw_recalc = mean(pval < 0.05),
                       size_centred_recalc = mean(pval_centered < 0.05),
                       n = .N), by = .(test, K)]
cat("\n== size sanity check (must match Table S7) ==\n")
print(size_check[order(K, test)])

# by construction the empirical-critical-value size is 5%; report the
# seed-matched null rejection under the pooled critical value as context
xcheck <- merge(size, crit_pool, by = c("test", "K"))[
  , .(null_rej_at_pooled_crit = mean(stat > crit_pool)), by = .(test, K)]
cat("\n== null rejection at pooled empirical crit (should be ~0.05) ==\n")
print(xcheck[order(K, test)])

if (HAS_POWER) {
  write_result(adj, paste0("e4_power_sizeadj", SFX), cfg)
  cat("\n== size-adjusted power at K = K_true ==\n")
  print(adj[K == K_true], nrows = 200L)
}

# ------------------------- (ii) contamination-refit by seed -------------------
if (HAS_POWER && !is.null(x$r2) && nrow(x$r2)) {
r2 <- x$r2[metric == "dev" & alt == "contamination" & K == K_true]
by_seed <- r2[, .(r2_micro = mean(r2_micro), r2_macro = mean(r2_macro),
                  n = .N), by = .(strength, train_seed)][order(strength,
                                                               train_seed)]
agg <- by_seed[, .(micro_mean = mean(r2_micro), micro_sd = sd(r2_micro),
                   macro_mean = mean(r2_macro), macro_sd = sd(r2_macro)),
               by = strength]
wide <- dcast(by_seed, train_seed ~ strength, value.var = "r2_micro")
write_result(by_seed, paste0("e4_contamrefit_by_seed", SFX), cfg)

cat("\n== contamination (refit): across-seed aggregate at K = K_true ==\n")
print(agg)
cat("\n== per-seed Micro R2 by strength (columns) ==\n")
print(wide)
mono <- wide[, mean(`0.05` > `0.1` & `0.1` < `0.2`)]
cat(sprintf(paste0(
  "\nREPORT: share of seeds with the non-monotone pattern ",
  "(dip at s=0.10, rebound at s=0.20): %.2f\n"), mono))
cat(paste0("If the dip-and-rebound holds across most seeds it is a systematic",
           " feature of the\nrefit alternative (the contamination topic",
           " becomes easy to fit at high strength),\nnot a run artifact;",
           " if driven by one or two seeds, it is Monte Carlo noise.\n"))
} # HAS_POWER && r2
log_msg("postprocess_e4 done for %s", basename(f))
