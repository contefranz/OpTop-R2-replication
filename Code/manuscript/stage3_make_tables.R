# Stage 3 (7 Oct 2026): LaTeX rows for the new supplement tables, generated from the Stage 3 CSVs
# (no model fitting). Inputs: Results/manuscript/stage3/stage3_*.csv and the cached E1 selection replicates.
suppressPackageStartupMessages({library(data.table); library(qs2)})
SA <- here::here("Results", "manuscript", "stage3"); OUT <- SA
tr <- fread(file.path(SA, "stage3_e2_reference_curves.csv"))
e2 <- fread(file.path(SA, "stage3_e2_falsecert.csv"))
e1 <- fread(file.path(SA, "stage3_e1_falsecert.csv"))
sr <- qs2::qs_read(here::here("Data/E1/e1_selreps_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2"))
sa <- as.data.table(sr$sel_all)[metric == "dev" & eval == "ho_reconstruction"]
# ---- Panel A: per fit ----
fits <- tr[, .(A40 = A[K == 40], A50 = A[K == 50], Kad01 = min(K[!is.na(A) & A <= 0.01]), Kad005 = min(K[!is.na(A) & A <= 0.005])), by = train_seed]
tg <- sa[rule == "total_gain" & eps == 0.01, .N, by = .(replicate, K_hat)][order(replicate, K_hat)]
fits <- merge(fits, tg[, .(sel = paste(sprintf("%d (%d)", K_hat, N), collapse = ", ")), by = .(train_seed = replicate)], by = "train_seed")
fwrite(fits, file.path(SA, "stage3_falsecert_fits.csv"))
rowsA <- fits[, sprintf("%d & %.4f & %.4f & %d & %d & %s \\\\", train_seed, A40, A50, Kad01, Kad005, sel)]
# ---- Panel B: counts per cell ----
rn <- c(total_gain = "Total gain, simultaneous", adjacent_simultaneous = "Adjacent gain, simultaneous", adjacent_pointwise = "Adjacent gain, pointwise")
b2 <- e2[, .(rule, eps, J_ev, certified, kstar = round(share_Kstar * 100), adequate = round(share_adequate * 100), false_cert)]
b1 <- e1[, .(rule, eps, J_ev = 5000L, certified, kstar = round(share_Kstar * 100), adequate = round(share_adequate * 100), false_cert)]
B <- rbind(b2, b1)[order(factor(rule, levels = names(rn)), -eps, J_ev)]
fwrite(B, file.path(SA, "stage3_falsecert_table.csv"))
rowsB <- B[, sprintf("%s & %.3f & %s & %d & %d & %d & %d \\\\", rn[rule], eps, format(J_ev, big.mark = ",", trim = TRUE), certified, kstar, adequate, false_cert)]
writeLines(c("% Panel A rows", rowsA, "% Panel B rows", rowsB), file.path(OUT, "tabS5_rows.tex"))
# ---- completion-residual table ----
ct <- fread(file.path(SA, "stage3_mdna_completion_tests.csv"))
tn <- c(T1_freq_contrast = "Test 1", T2_freq_strata = "Test 2", T3_fit_strata = "Test 3")
fp <- function(p) ifelse(p >= 1e-3, formatC(signif(p, 2), format = "fg", flag = "#"), { e <- floor(log10(p)); sprintf("%.1f\\times10^{%d}", p / 10^e, e) })
ct <- ct[order(K, factor(test, levels = names(tn)), factor(protocol, levels = c("reconstruction", "completion")))]
rowsC <- ct[, sprintf("%s & %d & %s & %.1f & $%s$ & %.1f & $%s$ \\\\", tn[test], K, protocol, W_iid, fp(p_iid), W_cluster, fp(p_cluster))]
writeLines(rowsC, file.path(OUT, "tabS25_rows.tex"))
ms <- fread(file.path(SA, "stage3_mdna_completion_masses.csv"))[test == "T2_freq_strata"]
cat("H (pp) by protocol/K:\n"); print(unique(ms[, .(protocol, K, H = round(H_pp, 3))]))
st <- fread(file.path(SA, "stage3_mdna_completion_strata.csv"))[test == "T2_freq_strata"]
cat("completion Test 2 clustered t:\n"); print(dcast(st, K + stratum ~ protocol, value.var = "t_cluster"), digits = 3)
# ---- present-firm rows ----
ss <- fread(file.path(SA, "stage3_mdna_subsample_selections.csv"))[sample == "present"]
g <- function(r, p, e) { v <- ss[rule == r & protocol == p & eps == e, K_hat]; if (length(v) != 1L || is.na(v)) "none" else as.character(v) }
rowsP <- sapply(names(rn), function(r) sprintf("Firms present in training, firm-clustered s.e. & %s & %s & %s & %s & %s \\\\", rn[r], g(r, "com", .01), g(r, "com", .005), g(r, "rec", .01), g(r, "rec", .005)))
writeLines(rowsP, file.path(OUT, "present_rows.tex"))
cat(readLines(file.path(OUT, "tabS5_rows.tex")), sep = "\n"); cat("\n"); cat(rowsC, sep = "\n"); cat("\n"); cat(rowsP, sep = "\n")
