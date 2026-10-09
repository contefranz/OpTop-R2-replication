# Stage 3 (7 Oct 2026): selections of the E2 coverage experiment judged against each training fit's
# reference curve (Table S5, Section 5.1). Reads the E2 rev2 object only; no fitting.
# Outputs: Results/manuscript/stage3/stage3_e2_{falsecert,reference_curves}.csv. Run first (stage3_e1_selreps.R reads its output).
suppressPackageStartupMessages({library(qs2); library(data.table)})
OUT <- here::here("Results", "manuscript", "stage3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
e2 <- qs2::qs_read(here::here("Data/E2/e2_results_E2_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2"))
tr <- as.data.table(e2$truth)[order(train_seed, K)]
# per-fit true remaining gain A_t(K) = max_{K'>K} mu_t(K') - mu_t(K); boundary K=100 undefined
tr[, A := { m <- mu; sapply(seq_along(m), function(i) if (i < length(m)) max(m[(i+1):length(m)]) - m[i] else NA_real_) }, by = train_seed]
tr[, A_se_rough := { s <- mu_se; sapply(seq_along(s), function(i) if (i < length(s)) sqrt(s[i]^2 + max(s[(i+1):length(s)])^2) else NA_real_) }, by = train_seed]
cat("Per-fit reference remaining gain A_t(K) (reconstruction, 20,000 reference docs per fit):\n")
print(dcast(tr[K >= 40], K ~ train_seed, value.var = "A"), digits = 3)
cat("\nMean reference curve mu(K) across fits and mean A(K):\n"); print(tr[, .(mu = mean(mu), A = mean(A), A_min = min(A), A_max = max(A)), by = K], digits = 4)
cat("\nSmallest K with A_t(K) <= eps, per fit:\n")
print(tr[, .(K_adequate_0.01 = min(K[!is.na(A) & A <= 0.01]), K_adequate_0.005 = min(K[!is.na(A) & A <= 0.005])), by = train_seed])
sa <- as.data.table(e2$sel_all)[metric == "dev"]
sa <- merge(sa, tr[, .(train_seed, K_hat = K, A_at_Khat = A)], by = c("train_seed", "K_hat"), all.x = TRUE)
res <- sa[, .(n = .N, certified = sum(certified), share_Kstar = mean(certified & K_hat == 40),
              share_adequate = mean(certified & A_at_Khat <= eps),
              false_cert = sum(certified & A_at_Khat > eps),
              false_cert_rate = mean(certified & A_at_Khat > eps),
              mean_A_at_Khat = mean(A_at_Khat[certified])), by = .(rule, eps, J_ev)][order(rule, eps, J_ev)]
cat("\nE2 selections vs per-fit reference curves (dev; 10 fits x 10 eval replications per cell):\n"); print(res, digits = 3)
# which false certifications, if any
fc <- sa[certified & A_at_Khat > eps, .(rule, eps, J_ev, train_seed, rep, K_hat, A_at_Khat)]
cat("\nFalse certifications (list):\n"); print(fc)
fwrite(res, file.path(OUT, "stage3_e2_falsecert.csv"))
fwrite(tr[, .(train_seed, K, mu, mu_se, A)], file.path(OUT, "stage3_e2_reference_curves.csv"))
