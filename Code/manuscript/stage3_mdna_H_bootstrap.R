# Stage 3 (7 Oct 2026): firm-cluster bootstrap for the displaced mass H at K = 50 and 180, from the
# centred firm sums behind main Figure 4 (Results/manuscript/residual_cluster_sums.csv). No model fitting.
rows <- read.csv(here::here("Results", "manuscript", "residual_cluster_sums.csv")); mass <- read.csv(here::here("Results", "manuscript", "residual_mass_data.csv"))
n <- 3466; G <- 3059; corr <- G/(G-1); B <- 20000; set.seed(20261007)
idx <- matrix(sample.int(G, G*B, replace = TRUE), nrow = B)
Hb <- function(K) { M <- as.matrix(rows[rows$K == K, c("m1","m2","m3","m4","m5")]); mk <- mass[mass$K == K, ]; mbar <- mk$mass_pp[order(mk$stratum)]
  list(H = 0.5*sum(abs(mbar)), draws = apply(idx, 1, function(ii) 0.5*sum(abs(mbar + colSums(M[ii, , drop = FALSE])/n)))) }
h50 <- Hb(50); h180 <- Hb(180); d <- h50$draws - h180$draws
out <- data.frame(quantity = c("H_K50", "H_K180", "H_K50_minus_H_K180"),
  estimate = c(h50$H, h180$H, h50$H - h180$H),
  boot_se = c(sd(h50$draws), sd(h180$draws), sd(d)),
  lwr95 = c(quantile(h50$draws, .025), quantile(h180$draws, .025), quantile(d, .025)),
  upr95 = c(quantile(h50$draws, .975), quantile(h180$draws, .975), quantile(d, .975)), B = B, seed = 20261007)
write.csv(out, here::here("Results", "manuscript", "stage3", "stage3_mdna_H_bootstrap.csv"), row.names = FALSE); print(out, digits = 4)
