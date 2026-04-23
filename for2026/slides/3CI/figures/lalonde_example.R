## LaLonde / Dehejia-Wahba running example for Stat 256 SOO_part2.
## Reproduces figures/lalonde_balance.pdf, lalonde_pscore.pdf,
## and figures/lalonde_estimates.pdf plus the estimates table
## figures/lalonde_estimates.rds.
##
## CI conventions (fixed-sample / model-based framing: nuisance fits
## treated as inputs, only Y random; SATT variance):
##   * DIM                              Neyman
##   * Stratified (u74 x u75)           cell-wise Neyman, ATT-aggregated
##   * Pooled OLS / Lin (ATE)           HC2 on lm()
##   * g-comp (ATE) via Lin             HC2 via the Lin identity
##                                      (arm-wise Neyman cross-check below)
##   * g-comp (ATT), any muhat_0 model  Neyman residual SE:
##                                        sd(Y - muhat_0 on treated) / sqrt(n_1)
##                                      (conditions on fitted muhat_0, same
##                                       "fix the nuisance" logic as the
##                                       weighted-LS-HC0 rows)
##   * Reweighting (PS-match, IPW,      weighted lm(Y ~ D + X, w = w),
##     ebal, prog-bal, kbal)              HC0 SE on the D coefficient
##
## Run from 3CI/ directory:  Rscript figures/lalonde_example.R

suppressPackageStartupMessages({
  library(kbal)
  library(ebal)
  library(Matching)
  library(sandwich)
  library(ggplot2)
  library(xgboost)
  library(dplyr)
  library(patchwork)
})

set.seed(20260421)

data(lalonde, package = "kbal")
covs  <- c("age", "educ", "black", "hisp", "married", "nodegr",
           "re74", "re75", "u74", "u75")
X     <- as.matrix(lalonde[, covs])
D     <- lalonde$nsw
Y     <- lalonde$re78
n     <- length(Y); n1 <- sum(D); n0 <- n - n1
dat   <- data.frame(Y = Y, D = D, X)
fm_pool    <- as.formula(paste("Y ~ D +", paste(covs, collapse = " + ")))
fm_outcome <- as.formula(paste("Y ~",      paste(covs, collapse = " + ")))

exp_bench <- 1794   # DW/LaLonde experimental-benchmark ATT, rounded.

## ---------- helpers -------------------------------------------------------

## Coefficient on D with HC0 SE, for any weighted or unweighted lm with X.
wls_D <- function(formula, data, w = NULL) {
  data$.w <- if (is.null(w)) rep(1, nrow(data)) else w
  m <- lm(formula, data = data, weights = .w)
  b <- coef(m)["D"]
  V <- sandwich::vcovHC(m, type = "HC0")
  c(est = unname(b), se = sqrt(V["D", "D"]))
}

## Coefficient on D with HC2 SE (used for pooled OLS and Lin).
ols_hc2_D <- function(formula, data) {
  m <- lm(formula, data = data)
  b <- coef(m)["D"]
  V <- sandwich::vcovHC(m, type = "HC2")
  c(est = unname(b), se = sqrt(V["D", "D"]))
}

## ---------- 1. DIM (Neyman) -----------------------------------------------
est_dim <- mean(Y[D == 1]) - mean(Y[D == 0])
se_dim  <- sqrt(var(Y[D == 1]) / n1 + var(Y[D == 0]) / n0)

## ---------- 2. Stratified on (u74, u75), ATT-aggregated --------------------
strat_key <- paste0(lalonde$u74, lalonde$u75)
cells     <- split(seq_len(n), strat_key)
cell_tab  <- do.call(rbind, lapply(cells, function(idx) {
  g1 <- idx[D[idx] == 1]; g0 <- idx[D[idx] == 0]
  if (length(g1) == 0 || length(g0) == 0) return(c(NA, NA, 0))
  dim_k <- mean(Y[g1]) - mean(Y[g0])
  var_k <- var(Y[g1]) / length(g1) + var(Y[g0]) / length(g0)
  w_k   <- length(g1) / n1
  c(dim_k, var_k, w_k)
}))
colnames(cell_tab) <- c("dim", "var", "w")
est_strat <- sum(cell_tab[, "w"] * cell_tab[, "dim"], na.rm = TRUE)
se_strat  <- sqrt(sum(cell_tab[, "w"]^2 * cell_tab[, "var"], na.rm = TRUE))

## ---------- 3. Pooled OLS (HC2) -------------------------------------------
h_pool    <- ols_hc2_D(fm_pool, dat)
est_pool  <- h_pool["est"]; se_pool <- h_pool["se"]

## ---------- 4. Lin: ATE and ATT (HC2) -------------------------------------
lin_est <- function(centre_vec) {
  Xc <- sweep(X, 2, centre_vec)
  colnames(Xc) <- paste0(covs, "_c")
  d  <- data.frame(Y = Y, D = D, Xc)
  rhs <- paste(colnames(Xc), collapse = " + ")
  fm  <- as.formula(paste("Y ~ D * (", rhs, ")"))
  ols_hc2_D(fm, d)
}
h_linATE <- lin_est(colMeans(X));           est_linATE <- h_linATE["est"]; se_linATE <- h_linATE["se"]
h_linATT <- lin_est(colMeans(X[D == 1, ])); est_linATT <- h_linATT["est"]; se_linATT <- h_linATT["se"]

## ---------- 5. g-computation (T-learner): ATE and ATT ---------------------
## Fit arm-wise OLS, impute, average.  With linear mu_d, the g-comp point
## estimate equals the D-coefficient of the Lin regression (see §4) when X
## is centred at the estimand's target mean; the HC2 SE on D from that Lin
## fit is therefore a valid SE for g-comp.  We also compute the equivalent
## arm-wise Neyman form as a cross-check.
m1 <- lm(fm_outcome, data = dat, subset = D == 1)
m0 <- lm(fm_outcome, data = dat, subset = D == 0)
mu1_all <- predict(m1, newdata = dat)
mu0_all <- predict(m0, newdata = dat)
est_gATE <- mean(mu1_all - mu0_all)
est_gATT <- mean((mu1_all - mu0_all)[D == 1])
se_gATE  <- se_linATE   # identity with Lin-ATE (X centred at overall mean)
## g-comp (ATT) SE: Neyman form on treated residuals, conditioning on
## fitted muhat_0.  Matches the "nuisance-as-input" convention used for
## the reweighting rows.
se_gATT  <- sd((Y - mu0_all)[D == 1]) / sqrt(n1)

## Cross-check: arm-wise Neyman variance.  mu_d_hat(x*) = (1, x*)' beta_d,
## so Var(mu_d_hat(x*)) = (1, x*)' V_d (1, x*) where V_d is arm-d HC2 vcov.
## Estimand-specific evaluation points:
##   ATE: x* = overall mean of X
##   ATT: x* = treated mean of X
arm_var <- function(armfit, target) {
  x_e <- c(1, as.numeric(target))
  V   <- sandwich::vcovHC(armfit, type = "HC2")
  as.numeric(t(x_e) %*% V %*% x_e)
}
se_gATE_armwise <- sqrt(arm_var(m1, colMeans(X)) +
                        arm_var(m0, colMeans(X)))
se_gATT_armwise <- sqrt(arm_var(m1, colMeans(X[D == 1, ])) +
                        arm_var(m0, colMeans(X[D == 1, ])))
cat(sprintf(
  "g-comp ATE SE: Lin HC2 = %.1f   arm-wise Neyman (HC2) = %.1f\n",
  se_gATE, se_gATE_armwise))
cat(sprintf(
  "g-comp ATT SE: Lin HC2 = %.1f   arm-wise Neyman (HC2) = %.1f\n",
  se_gATT, se_gATT_armwise))

## ---------- 6. PS 1-to-1 NN matching (ATT), then weighted lm on (D, X) ----
ps_fit <- glm(D ~ ., data = data.frame(D = D, X), family = binomial)
ps     <- predict(ps_fit, type = "response")
mm     <- Match(Y = Y, Tr = D, X = ps, estimand = "ATT", M = 1, replace = TRUE)

w_match <- numeric(n)
w_match[D == 1] <- 1
for (k in seq_along(mm$index.control)) {
  j <- mm$index.control[k]
  w_match[j] <- w_match[j] + mm$weights[k]
}
h_psm <- wls_D(fm_pool, dat, w = w_match)
est_psm <- h_psm["est"]; se_psm <- h_psm["se"]

## ---------- 7. IPW (stabilized / HT form, ATT), weighted lm on (D, X) ----
## ATT weights (stabilized-scale): treated are 1, controls reweighted by
## pi/(1 - pi) to match the treated-X distribution.  A stabilization
## constant P(D=1)/P(D=0) on controls cancels under weighted LS, so we use
## the clean ATT form directly.
w_ipw <- ifelse(D == 1, 1, ps / (1 - ps))
h_ipw <- wls_D(fm_pool, dat, w = w_ipw)
est_ipw <- h_ipw["est"]; se_ipw <- h_ipw["se"]

## ---------- 8. Entropy balancing (ATT), weighted lm on (D, X) -------------
eb <- ebalance(Treatment = D, X = X, print.level = -1)
w_ebal <- numeric(n); w_ebal[D == 1] <- 1; w_ebal[D == 0] <- eb$w
h_ebal <- wls_D(fm_pool, dat, w = w_ebal)
est_ebal <- h_ebal["est"]; se_ebal <- h_ebal["se"]

## ---------- 9. Prognostic balancing (ATT) --------------------------------
## Workflow:
##   (1) fit muhat_0(X) on controls only;
##   (2) project muhat_0 onto all units (including treated);
##   (3) find control weights via ebalance so control mean of muhat_0
##       equals treated mean of muhat_0 (mean balance on the prog. score);
##   (4) weighted lm(Y ~ D + X, w), HC0 SE on D.
##
## Two flavours of the muhat_0 model: linear (lm) and flexible (xgboost).

## (a) lm prognostic model -------------------------------------------------
prog_lm <- mu0_all    # m0 above is an OLS fit on controls
eb_lm <- ebalance(Treatment = D, X = as.matrix(prog_lm), print.level = -1)
w_prog_lm <- numeric(n); w_prog_lm[D == 1] <- 1; w_prog_lm[D == 0] <- eb_lm$w
h_prog_lm <- wls_D(fm_pool, dat, w = w_prog_lm)
est_prog_lm <- h_prog_lm["est"]; se_prog_lm <- h_prog_lm["se"]

## (b) xgboost prognostic model -------------------------------------------
set.seed(20260421)
X_mat <- as.matrix(X)
dctrl <- xgb.DMatrix(data = X_mat[D == 0, ], label = Y[D == 0])
## Modest regularisation so the prog. score isn't overfit on 2,490 controls.
xgb_prog <- xgb.train(
  params = list(objective      = "reg:squarederror",
                eta            = 0.05,
                max_depth      = 4,
                min_child_weight = 10,
                subsample      = 0.8,
                colsample_bytree = 0.8),
  data   = dctrl,
  nrounds = 400,
  verbose = 0
)
prog_xgb <- predict(xgb_prog, xgb.DMatrix(data = X_mat))
eb_xgb <- ebalance(Treatment = D, X = as.matrix(prog_xgb), print.level = -1)
w_prog_xgb <- numeric(n); w_prog_xgb[D == 1] <- 1; w_prog_xgb[D == 0] <- eb_xgb$w
h_prog_xgb <- wls_D(fm_pool, dat, w = w_prog_xgb)
est_prog_xgb <- h_prog_xgb["est"]; se_prog_xgb <- h_prog_xgb["se"]

## ---------- 9b. g-comp (ATT) with xgboost mu_0 ---------------------------
## Form 1 g-comp-ATT: tau = mean over treated of [Y_i - muhat_0_xgb(X_i)].
## Same Neyman residual SE as the lm g-comp (ATT), just a different muhat_0:
##   SE = sd(Y_i - muhat_0(X_i) on treated) / sqrt(n_1).
## Conditions on the fitted xgb muhat_0 (nuisance-as-input convention).
est_gATT_xgb <- mean((Y - prog_xgb)[D == 1])
se_gATT_xgb  <- sd((Y - prog_xgb)[D == 1]) / sqrt(n1)

## ---------- 10. Kernel balancing (kbal), ATT via weighted lm --------------
## Spec from Example 1 in ?kbal (no b set: use default search for b).
cat("Running kbal ... (may take a minute)\n")
xvars_kbal <- c("age","black","educ","hisp","married",
                "re74","re75","nodegr","u74","u75")
kb <- kbal(allx      = lalonde[, xvars_kbal],
           treatment = lalonde$nsw,
           fullSVD   = TRUE,
           printprogress = FALSE)
w_kbal <- kb$w
## kbal returns sample-wide weights; treated units already at 1/n1 scale.
h_kbal <- wls_D(fm_pool, dat, w = w_kbal)
est_kbal <- h_kbal["est"]; se_kbal <- h_kbal["se"]

## ---------- collect + save ------------------------------------------------
## Labels: reserve parentheses for the estimand (ATE / ATT).  Methods that
## don't cleanly target either (naive DIM, pooled OLS under heterogeneity)
## are left unparenthesised.
res <- data.frame(
  method = c("DIM", "Stratify u74,u75 (ATT)",
             "Pooled OLS",
             "g-comp / Lin (ATE)",
             "g-comp lm (ATT)",
             "g-comp xgb (ATT)",
             "PS match 1:1 (ATT)",
             "IPW stab. (ATT)",
             "Entropy balancing (ATT)",
             "Prognostic bal lm (ATT)",
             "Prognostic bal xgb (ATT)",
             "kbal (ATT)"),
  estimate = c(est_dim, est_strat,
               est_pool,
               est_gATE,
               est_gATT,
               est_gATT_xgb,
               est_psm, est_ipw,
               est_ebal,
               est_prog_lm, est_prog_xgb,
               est_kbal),
  se = c(se_dim, se_strat,
         se_pool,
         se_gATE,
         se_gATT,
         se_gATT_xgb,
         se_psm, se_ipw,
         se_ebal,
         se_prog_lm, se_prog_xgb,
         se_kbal)
)
res$lo <- res$estimate - 1.96 * res$se
res$hi <- res$estimate + 1.96 * res$se
## Sort worst (most negative point estimate) at top of plot to best at bottom.
res <- res[order(res$estimate), ]
res$method <- factor(res$method, levels = rev(res$method))

cat("\nEstimates (point, SE, 95% CI):\n")
print(res, digits = 4, row.names = FALSE)

saveRDS(res, "figures/lalonde_estimates.rds")

## ---------- Figure A: estimates comparison --------------------------------
## Heavy x = 0 line; lighter dashed benchmark line.
p_est <- ggplot(res, aes(y = method, x = estimate)) +
  geom_vline(xintercept = 0,        linetype = "solid",
             colour = "grey20", linewidth = 0.9) +
  geom_vline(xintercept = exp_bench, linetype = "dashed",
             colour = "#64b5f6", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, na.rm = TRUE,
                 colour = "grey30") +
  geom_point(size = 2.4, colour = "black") +
  scale_x_continuous(labels = function(x)
    paste0("$", formatC(x, format = "d", big.mark = ","))) +
  labs(x = "Estimated effect on 1978 earnings (USD)", y = NULL,
       title = "LaLonde / Dehejia-Wahba: observational estimates",
       subtitle = sprintf("NSW treated vs PSID controls.  Experimental benchmark ~ $%d (light blue dashed).", exp_bench)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave("figures/lalonde_estimates.pdf", p_est, width = 7.4, height = 5.1)

## ---------- Figure B: pscore overlap --------------------------------------
## Two x-panels.  Left zooms on Pr(D=1|X) < 0.05 where 91% of controls live
## (only 9 treated there, so ATE-style reweighting would make each of those
## 9 carry the weight of hundreds of controls).  Right shows Pr(D=1|X) >=
## 0.05 where counts are within an order of magnitude and ATT is workable.
## We use explicit count tables (not geom_histogram) so we can annotate
## individual bins (treated counts as a blue strip at top of the zoom panel,
## since the actual treated bars are too small to see next to the 1,854-
## tall control bar in the first bin).

mk_counts <- function(ps, D, brk) {
  keep <- ps >= min(brk) & ps <= max(brk)
  ps1 <- ps[keep]; D1 <- D[keep]
  bin_idx <- findInterval(ps1, brk, rightmost.closed = TRUE,
                          all.inside = TRUE)
  tab <- as.data.frame(table(
    bin   = factor(bin_idx, levels = seq_len(length(brk) - 1)),
    group = factor(D1, levels = c(1, 0),
                   labels = c("NSW treated", "PSID control"))
  ), stringsAsFactors = FALSE)
  tab$bin <- as.integer(tab$bin)
  tab$lo  <- brk[tab$bin]
  tab$hi  <- brk[tab$bin + 1]
  tab$mid <- (tab$lo + tab$hi) / 2
  tab$group <- factor(tab$group,
                      levels = c("NSW treated", "PSID control"))
  names(tab)[names(tab) == "Freq"] <- "n"
  tab
}

brk_zoom <- seq(0, 0.05, by = 0.005)    # 10 bins of width 0.005
brk_rest <- seq(0.05, 1.00, by = 0.025) # 38 bins of width 0.025

cz <- mk_counts(ps, D, brk_zoom)
cr <- mk_counts(ps, D, brk_rest)

col_trt <- "#1565c0"; col_ctl <- "#c62828"
pal <- c("NSW treated" = col_trt, "PSID control" = col_ctl)

cz_ctl_lab <- cz |> dplyr::filter(group == "PSID control", n > 0) |>
  dplyr::mutate(label = format(n, big.mark = ","))
cz_trt_lab <- cz |> dplyr::filter(group == "NSW treated") |>
  dplyr::mutate(label = as.character(n))

n_tot_ctl_zoom <- sum(D == 0 & ps <  0.05)
n_tot_trt_zoom <- sum(D == 1 & ps <  0.05)
n_tot_ctl_rest <- sum(D == 0 & ps >= 0.05)
n_tot_trt_rest <- sum(D == 1 & ps >= 0.05)

y_max_zoom <- max(cz$n)

p_zoom <- ggplot(cz, aes(x = mid, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.0045), width = 0.004,
           colour = "grey20", linewidth = 0.15) +
  geom_text(data = cz_ctl_lab,
            aes(x = mid, y = n, label = label),
            vjust = -0.4, size = 2.6, fontface = "bold",
            colour = col_ctl, inherit.aes = FALSE) +
  geom_text(data = cz_trt_lab,
            aes(x = mid, label = label),
            y = y_max_zoom * 1.18, size = 3.0, fontface = "bold",
            colour = col_trt, inherit.aes = FALSE) +
  annotate("text", x = 0.025, y = y_max_zoom * 1.28,
           label = "treated count per bin (bars too small to see):",
           size = 2.8, fontface = "italic", colour = col_trt) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(-0.002, 0.052),
                     breaks = seq(0, 0.05, 0.01),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "Pr(D=1 | X)",
       y = "Count per 0.005 bin",
       title = sprintf("Pr(D=1|X) < 0.05:  %d treated for %s controls",
                       n_tot_trt_zoom,
                       format(n_tot_ctl_zoom, big.mark = ",")),
       fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

p_rest <- ggplot(cr, aes(x = mid, y = n, fill = group)) +
  geom_col(position = position_dodge(width = 0.022), width = 0.02,
           colour = "grey20", linewidth = 0.15) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(limits = c(0.05, 1.0), breaks = seq(0.1, 1.0, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Pr(D=1 | X)",
       y = "Count per 0.025 bin",
       title = sprintf("Pr(D=1|X) >= 0.05:  %d treated for %d controls",
                       n_tot_trt_rest, n_tot_ctl_rest),
       fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = c(0.5, 0.97), legend.justification = c(0.5, 1),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", colour = NA),
        plot.title = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

p_ps <- p_zoom + p_rest + plot_layout(widths = c(1, 2.2)) +
  plot_annotation(
    title = "P-score overlap: NSW treated vs PSID controls",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))

ggsave("figures/lalonde_pscore.pdf", p_ps, width = 10.2, height = 4.2)

## ---------- Figure C: pre / post balance (entropy balancing) --------------
smd <- function(x, D, w = NULL) {
  if (is.null(w)) w <- rep(1, length(D))
  m1 <- sum(w[D == 1] * x[D == 1]) / sum(w[D == 1])
  m0 <- sum(w[D == 0] * x[D == 0]) / sum(w[D == 0])
  s1 <- sqrt(sum(w[D == 1] * (x[D == 1] - m1)^2) / sum(w[D == 1]))
  (m1 - m0) / s1
}
smd_pre  <- sapply(covs, function(v) smd(X[, v], D))
smd_ipw  <- sapply(covs, function(v) smd(X[, v], D, w_ipw))
smd_ebal <- sapply(covs, function(v) smd(X[, v], D, w_ebal))
bal_df <- data.frame(
  covariate = factor(rep(covs, 3), levels = rev(covs)),
  smd       = c(smd_pre, smd_ipw, smd_ebal),
  when      = factor(rep(c("Unweighted",
                           "IPW (ATT)",
                           "Entropy balancing (ATT)"),
                         each = length(covs)),
                     levels = c("Unweighted",
                                "IPW (ATT)",
                                "Entropy balancing (ATT)")))
p_bal <- ggplot(bal_df, aes(x = smd, y = covariate,
                            colour = when, shape = when)) +
  geom_vline(xintercept = 0,        linetype = "dotted", colour = "grey50") +
  geom_vline(xintercept = c(-.1, .1), linetype = "dashed", colour = "grey80") +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("Unweighted"               = "#c62828",
                                 "IPW (ATT)"                = "#1565c0",
                                 "Entropy balancing (ATT)"  = "#2e7d32")) +
  labs(x = "Standardised mean difference (treated - control)",
       y = NULL, colour = NULL, shape = NULL,
       title = "Covariate balance on LaLonde") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave("figures/lalonde_balance.pdf", p_bal, width = 6.8, height = 4.6)

## ---------- Figure D: curse of dimensionality --------------------------
## For X uniform in R^P with unit density, Euclidean radius r for a ball
## around a point to contain fraction f of the mass:
##   V_P * r^P = f   where   V_P = pi^(P/2) / Gamma(P/2 + 1)
## so  r = (f / V_P)^(1/P).  Teaches how fast "close" gets far as P grows.
f_capture <- 0.01
V_P <- function(P) pi^(P/2) / gamma(P/2 + 1)
P_vec <- c(1, 2, 4, 8, 16, 32, 64, 128)
r_by_P <- sapply(P_vec, function(P) (f_capture / V_P(P))^(1/P))
cod_df <- data.frame(P = P_vec, r = r_by_P)
cat("\nCurse-of-dim: Euclidean r to capture 1% of uniform mass vs P\n")
print(cod_df, row.names = FALSE)

p_cod <- ggplot(cod_df, aes(y = factor(P, levels = rev(P_vec)), x = r)) +
  geom_col(fill = "#1565c0", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", r)),
            hjust = -0.12, size = 3.2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(y = "dim(X) = P",
       x = "Euclidean distance to the 1% nearest neighbours",
       title = "Curse of dimensionality",
       subtitle = "If each X ~ unif(0,1), for dim(X) = P, how far away are the 1% nearest neighbours?") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave("figures/curse_of_dim.pdf", p_cod, width = 7.0, height = 4.0)

## ---------- Figure E: ATE-IPW weight distribution, treated arm ----------
## Stabilised ATE weights: w_i = P(D=1)/pi for treated, P(D=0)/(1-pi) for controls.
## Show treated arm only (controls' ESS is fine: many discarded but plenty left).
pD1 <- mean(D); pD0 <- 1 - pD1
w_ate <- ifelse(D == 1, pD1 / ps, pD0 / (1 - ps))
ess_arm <- function(w) (sum(w))^2 / sum(w^2)
ess_T <- ess_arm(w_ate[D == 1])
ipw_df <- data.frame(w = w_ate[D == 1])
p_ipw <- ggplot(ipw_df, aes(x = w)) +
  geom_histogram(bins = 30, fill = "#1565c0", colour = "white",
                 linewidth = 0.2) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(x = "Stabilised ATE-IPW weight  (log scale)", y = "Count",
       title = sprintf("ATE-IPW weights on LaLonde, NSW treated  (n = %d, ESS = %.0f)",
                       n1, ess_T)) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave("figures/lalonde_ipw_weights.pdf", p_ipw, width = 6.8, height = 3.4)

cat("\nDone. Figures in figures/.\n")
