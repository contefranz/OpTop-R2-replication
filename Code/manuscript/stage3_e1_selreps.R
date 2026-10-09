# Stage 3 (7 Oct 2026): E1 selection replicates judged against the per-fit reference curves (Table S5,
# J_ev = 5,000 rows) and across-fit Monte Carlo SEs of the Table 2 shares. Needs stage3_e2_falsecert.R output.
# Outputs: Results/manuscript/stage3/stage3_e1_{table2_mcse,falsecert}.csv.
suppressPackageStartupMessages({library(qs2); library(data.table)})
SP <- here::here("Results", "manuscript", "stage3")
tr <- fread(file.path(SP, "stage3_e2_reference_curves.csv"))
sr <- qs2::qs_read(here::here("Data/E1/e1_selreps_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2"))
cat("selreps names:", names(sr), "\n")
sa <- as.data.table(sr$sel_all)[metric == "dev"]
cat("sel_all cols:", names(sa), "; rows:", nrow(sa), "; evals:", unique(sa$eval), "; replicates:", length(unique(sa$replicate)), "; reps per replicate:", length(unique(sa$rep)), "\n")
# sanity: E1 fit identity with E2 train_seed -> compare reconstruction Macro at K=40 per replicate (first eval corpus) with E2 truth mu
su <- as.data.table(sr$summary)[metric == "dev" & K == 40]
chk <- merge(su[, .(r2 = mean(r2_macro)), by = .(eval, replicate)], tr[K == 40, .(replicate = train_seed, mu_ref = mu)], by = "replicate")
cat("\nFit identity check (K=40): E1 selreps mean Macro per replicate vs E2 reference mu (reconstruction):\n"); print(dcast(chk, replicate + mu_ref ~ eval, value.var = "r2"), digits = 4)
# across-fit Monte Carlo SE of the Table 2 shares
tab <- sa[, .(exact = mean(certified & K_hat == 40)), by = .(eval, rule, eps, replicate)][, .(share = mean(exact), mc_se_fits = sd(exact)/sqrt(.N), n_fits = .N), by = .(eval, rule, eps)][order(eval, rule, eps)]
cat("\nTable 2 shares with across-fit Monte Carlo SE (10 fits x 10 corpora):\n"); print(tab, digits = 3)
# false certification against the per-fit E2 reference curve (reconstruction only)
rec <- sa[eval %in% grep("rec", unique(sa$eval), value = TRUE)]
rec <- merge(rec, tr[, .(replicate = train_seed, K_hat = K, A = A)], by = c("replicate", "K_hat"), all.x = TRUE)
res <- rec[, .(n = .N, certified = sum(certified), share_Kstar = mean(certified & K_hat == 40), share_adequate = mean(certified & A <= eps), false_cert = sum(certified & A > eps), max_A_certified = max(A[certified], na.rm = TRUE)), by = .(rule, eps)][order(rule, eps)]
cat("\nE1 reconstruction selections (n=100) vs per-fit reference curves:\n"); print(res, digits = 3)
# per-fit view for total gain eps=0.01: which K_hat, and A(K*=40) per fit
pf <- rec[rule == "total_gain" & eps == 0.01, .(K_hats = paste(sort(unique(K_hat)), collapse = "/"), n_cert = sum(certified)), by = replicate]
pf <- merge(pf, tr[K == 40, .(replicate = train_seed, A40 = round(A, 4))], by = "replicate")
cat("\nTotal gain eps=0.01 reconstruction: selections by fit and the fit's A(40):\n"); print(pf)
# comparators in the selreps object?
if (!is.null(sr$khat)) { kh <- as.data.table(sr$khat); cat("\nkhat rules:", unique(kh$rule), "\n") }
fwrite(tab, file.path(SP, "stage3_e1_table2_mcse.csv")); fwrite(res, file.path(SP, "stage3_e1_falsecert.csv"))
