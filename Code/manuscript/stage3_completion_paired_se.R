# Stage 3 (9 Oct 2026): paired standard errors of the per-fit COMPLETION remaining gain A_D(K) for
# K = 40 and 50 (Supplement Section S5.1 and Table S5 notes: "In completion fit 7, the paired standard
# error for the comparison attaining the estimated A_D(40) = 0.00474 is about 0.00086"). For each
# training fit, the attaining comparison is K -> argmax_{K' > K} mu(K'); its paired gain and standard
# error come from the project's paired_gains_all() on the stored reference-document scores.
# Input : Data/E2/e2_completion_truth_rev2.qs2 (written by stage3_completion_reference.R full;
#         Zenodo results record). No fitting, no rescoring.
# Output: Results/manuscript/stage3/stage3_completion_paired_se.csv
# Usage : Rscript Code/manuscript/stage3_completion_paired_se.R      (from the package root; ~10 s)
suppressPackageStartupMessages(source(here::here("Code", "R", "source_all.R")))
OUT <- here::here("Results", "manuscript", "stage3"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
tr <- qs_read(p_data("E2", "e2_completion_truth_rev2.qs2"))
stopifnot(identical(tr$mode, "full"))
cur <- as.data.table(tr$curves); d_all <- as.data.table(tr$doc)
res <- rbindlist(lapply(sort(unique(cur$train_seed)), function(t) {
  ct <- cur[train_seed == t]
  pg <- paired_gains_all(d_all[train_seed == t & metric == "dev"], 0.05, NULL, "all_pairs")
  rbindlist(lapply(c(40L, 50L), function(k0) {
    k_att <- ct[K > k0][which.max(mu), K]
    row <- pg[K == k0 & K_to == k_att]
    stopifnot(nrow(row) == 1L, abs(row$delta_mean - ct[K == k0, A]) < 1e-12)   # paired gain = A_D(k0)
    data.table(train_seed = t, K = k0, K_attaining = k_att, A_D = ct[K == k0, A],
               paired_se = row$se, n_docs = row$n)
  }))
}))
print(res, digits = 4)
f7 <- res[train_seed == 7L & K == 40L]
cat(sprintf("\nFit 7, A_D(40) = %.5f via 40 -> %d; paired SE = %.5f (manuscript: about 0.00086)\n",
            f7$A_D, f7$K_attaining, f7$paired_se))
fwrite(res, file.path(OUT, "stage3_completion_paired_se.csv"))
