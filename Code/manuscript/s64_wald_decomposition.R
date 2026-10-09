# Stage 3 (7 Oct 2026): decomposition of the firm-clustered Test 2 Wald statistic of the MD&A
# application at K = 50 and 180 (Supplement Section S6.4; main Section 6.3: "equal-weight mean
# contrast t = 3.2", homogeneity "p = 0.08"). Works from the anonymous centred firm sums behind main
# Figure 4 (written by make_main_figures.R); no model fitting, no rescoring.
#   * joint clustered Wald statistic of the four contrasts (f_b vs f_5) and the sum of squared marginal t
#   * joint test of the three lowest groups; homogeneity of the four lower groups
#   * GLS common-shift statistic, with W = W_homogeneity + W_common_shift exactly
#   * all pairwise and equal-weight mean contrasts; Scheffe bound
# Inputs : Results/manuscript/{residual_cluster_sums,residual_contrast_data,residual_mass_data}.csv
# Output : Results/manuscript/stage3/stage3_mdna_s64_wald_subtests.csv (and the printed report)
# Usage  : Rscript Code/manuscript/s64_wald_decomposition.R      (from the package root; < 1 s)
suppressPackageStartupMessages(library(data.table))
IN <- here::here("Results", "manuscript"); OUT <- file.path(IN, "stage3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rows <- read.csv(file.path(IN, "residual_cluster_sums.csv"))
con  <- read.csv(file.path(IN, "residual_contrast_data.csv"))
mass <- read.csv(file.path(IN, "residual_mass_data.csv"))
n <- 3466; G <- 3059; corr <- G/(G-1)
res <- list()
add <- function(K, quantity, value, df = NA_integer_, p = NA_real_)
  res[[length(res) + 1L]] <<- data.table(K = K, quantity = quantity, value = value, df = df, p_value = p)
for (K in c(50,180)) {
  C <- as.matrix(rows[rows$K==K, c('c1','c2','c3','c4')])
  M <- as.matrix(rows[rows$K==K, c('m1','m2','m3','m4','m5')])
  cat(sprintf("\n=== K=%d: clusters=%d; colSums c: %s ; colSums m: %s\n", K, nrow(C),
      paste(signif(colSums(C),3),collapse=' '), paste(signif(colSums(M),3),collapse=' ')))
  m <- con$gbar[con$K==K]
  V <- corr*crossprod(C)/n^2
  W <- drop(t(m)%*%solve(V,m)); cat("Test2 clustered W =",round(W,2)," p=",signif(pchisq(W,4,lower=FALSE),2),"\n")
  add(K, "W_test2_cluster", W, 4L, pchisq(W, 4, lower = FALSE))
  se <- sqrt(diag(V)); tt <- m/se; cat("clustered t:",round(tt,2)," sum t^2 =",round(sum(tt^2),1),"\n")
  for (b in 1:4) add(K, sprintf("t_f%d_vs_f5", b), tt[b])
  add(K, "sum_t_squared", sum(tt^2))
  R <- V/outer(se,se); cat("correlations:",round(sort(R[upper.tri(R)]),3),"\n")
  add(K, "min_contrast_correlation", min(R[upper.tri(R)])); add(K, "max_contrast_correlation", max(R[upper.tri(R)]))
  W3 <- drop(t(m[1:3])%*%solve(V[1:3,1:3],m[1:3])); cat("3 lowest joint W =",round(W3,2)," p=",round(pchisq(W3,3,lower=FALSE),3),"\n")
  add(K, "W_three_lowest", W3, 3L, pchisq(W3, 3, lower = FALSE))
  D <- rbind(c(1,-1,0,0),c(1,0,-1,0),c(1,0,0,-1)); d <- D%*%m
  WH <- drop(t(d)%*%solve(D%*%V%*%t(D),d)); cat("homogeneity of 4 lower strata W =",round(WH,2)," p=",signif(pchisq(WH,3,lower=FALSE),2),"\n")
  add(K, "W_homogeneity_lower4", WH, 3L, pchisq(WH, 3, lower = FALSE))
  nm <- paste0('f',1:5)
  cv <- function(a,b){c<-rep(0,4); if(a<5) c[a]<-c[a]+1; if(b<5) c[b]<-c[b]-1; c}
  for (a in 1:4) for (b in (a+1):5) { c<-cv(a,b); tab <- drop(c%*%m)/sqrt(drop(t(c)%*%V%*%c))
    cat(sprintf("  %s-%s t=%.2f\n", nm[a], nm[b], tab)); add(K, sprintf("t_%s_minus_%s", nm[a], nm[b]), tab) }
  one <- rep(1/4,4); tew <- drop(one%*%m)/sqrt(drop(t(one)%*%V%*%one)); cat("equal-weight mean t =",round(tew,2),"\n")
  add(K, "t_equal_weight_mean", tew)
  Vi <- solve(V); ones <- rep(1,4); Wg <- drop(t(ones)%*%Vi%*%m)^2/drop(t(ones)%*%Vi%*%ones)
  cat("GLS common-shift W =",round(Wg,2),"; W_H + W_gls =",round(WH+Wg,2)," vs W =",round(W,2),"; GLS shift t =",round(sqrt(Wg)*sign(drop(t(ones)%*%Vi%*%m)),2),"\n")
  stopifnot(abs(WH + Wg - W) < 1e-8 * W)   # exact additive split
  add(K, "W_gls_common_shift", Wg, 1L, pchisq(Wg, 1, lower = FALSE))
  add(K, "t_gls_common_shift", sqrt(Wg)*sign(drop(t(ones)%*%Vi%*%m)))
  cat("Scheffe bound sqrt(chi2_4,.95) =",round(sqrt(qchisq(.95,4)),3),"\n"); add(K, "scheffe_bound", sqrt(qchisq(.95, 4)))
  mk <- mass[mass$K==K,]; cat(sprintf("  group %d: mass %+.4f pp, se %.4f, per-word SD %.2e\n", mk$stratum, mk$mass_pp, mk$se, mk$se/100/mk$words), sep='')
  # cross-check: per-word group SDs from cluster sums of masses (centred, in pp)
  Vm <- corr*crossprod(M)/n^2; cat("  per-word SD from m-cols:", signif(sqrt(diag(Vm))/100/mk$words,3), "\n")
}
fwrite(rbindlist(res), file.path(OUT, "stage3_mdna_s64_wald_subtests.csv"))
