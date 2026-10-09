# Stage 3 (7 Oct 2026): MD&A selections and levels on the present-firm / absent-firm subsamples with
# firm-clustered standard errors (Table S28 rows), and firm-clustered Micro-Macro gap intervals (S6.2, Section 6.3).
# Cache-only. Outputs: Results/manuscript/stage3/stage3_mdna_{subsample_selections,subsample_levels,gap_cluster}.csv.
suppressPackageStartupMessages(source(here::here("Code", "R", "source_all.R")))
SP <- here::here("Results", "manuscript", "stage3"); EPS <- c(0.01, 0.005); ALPHA <- 0.05
x <- qs_read(p_data("MDNA", "mdna_results_MDNA_2015_2016_rev2.qs2"))
prep <- qs_read(p_data("MDNA", "mdna_prep_2015_2016.qs2"))
stopifnot(identical(prep$dv_ev$doc_id, rownames(prep$dtm_ev)))
cik <- data.table(doc_id = prep$dv_ev$doc_id, cluster_id = as.character(prep$dv_ev$cik))
in_tr <- prep$dv_ev$cik %in% prep$dv_train$cik
ids <- list(full = prep$dv_ev$doc_id, present = prep$dv_ev$doc_id[in_tr], absent = prep$dv_ev$doc_id[!in_tr])
cat("docs:", sapply(ids, length), "; firms:", sapply(ids, function(i) uniqueN(cik[doc_id %in% i, cluster_id])), "\n")
docs <- list(com = x$doc_com, rec = x$doc_rec)
rules_for <- function(doc_list, cluster = NULL) rbindlist(lapply(names(doc_list), function(pn) select_k_all_rules(doc_list[[pn]], EPS, ALPHA, cluster)[, protocol := pn]))
sel <- rbindlist(lapply(names(ids), function(s) rules_for(lapply(docs, function(d) d[doc_id %in% ids[[s]]]), cik)[, sample := s]))
sel <- sel[metric == "dev"]
cat("\nSelections (Deviance, firm-clustered s.e.):\n")
print(dcast(sel, sample + rule ~ protocol + eps, value.var = "K_hat")[order(factor(sample, levels = names(ids)), rule)])
# retained documents and index levels at K=50 and 180
lev <- rbindlist(lapply(names(ids), function(s) rbindlist(lapply(names(docs), function(pn) {
  d <- docs[[pn]][metric == "dev" & doc_id %in% ids[[s]]]
  ci <- macro_ci(d)[K %in% c(50, 180, 200), .(K, r2_macro, se, n)]
  allp <- paired_gains_all(d, ALPHA, cik, "all_pairs")
  ci[, `:=`(sample = s, protocol = pn, gain_to_200 = allp[K_to == 200][match(ci$K, K), delta_mean],
            gain_50_to_200_ub_simul = allp[K == 50 & K_to == 200, ub_simul])]
}))))
cat("\nMacro Deviance levels, retained n, estimated gain to K=200 (firm-clustered simultaneous UB for 50->200):\n"); print(lev[order(protocol, factor(sample, levels = names(ids)), K)], digits = 4)
# firm-clustered Micro-Macro gap intervals (delta-method influence function, clustered by firm)
gap_cl <- function(doc_dt, cl) {
  z <- qnorm(0.975)
  doc_dt[!is.na(r2_doc) & d_null > 0, {
    u <- r2_doc; v <- d_null; n <- .N; ub <- mean(u); vb <- mean(v); suv <- mean((u - ub) * (v - vb)); gap <- suv / vb
    psi <- ((u - ub) * (v - vb) - suv) / vb - gap * (v - vb) / vb
    g <- cl$cluster_id[match(doc_id, cl$doc_id)]; G <- uniqueN(g)
    s_iid <- sd(psi) / sqrt(n); s_cl <- sqrt(G / (G - 1) * sum(tapply(psi, g, sum)^2)) / n
    .(gap = gap, se_iid = s_iid, se_cluster = s_cl, ratio = s_cl / s_iid, n = n, G = G, lwr_cl = gap - z * s_cl, upr_cl = gap + z * s_cl)
  }, by = .(K, metric)]
}
gp <- rbindlist(lapply(names(docs), function(pn) gap_cl(docs[[pn]][metric == "dev"], cik)[, protocol := pn]))
cat("\nStored gap table eval labels:", unique(x$gap$eval), "\n")
chk <- merge(gp, as.data.table(x$gap)[metric == "dev", .(K, protocol = ifelse(grepl("com", eval), "com", "rec"), se_stored = se, gap_stored = gap)], by = c("K", "protocol"))
cat("Firm-clustered gap SE vs stored document-level SE (should match se_iid):\n"); print(chk[K %in% c(10, 50, 100, 180, 200), .(protocol, K, gap, gap_stored, se_iid, se_stored, se_cluster, ratio)], digits = 4)
cat("\nRatio cluster/iid SE for the gap over all K: range", range(gp$ratio), "\n")
fwrite(sel, file.path(SP, "stage3_mdna_subsample_selections.csv")); fwrite(lev, file.path(SP, "stage3_mdna_subsample_levels.csv")); fwrite(gp, file.path(SP, "stage3_mdna_gap_cluster.csv"))
