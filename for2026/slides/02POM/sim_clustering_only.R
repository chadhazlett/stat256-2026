# sim_clustering_only.R
# Produces four figures for the Stat 256 experiments lecture, showing the
# sampling distribution of tauhat under four designs/estimators:
#   sim_simple.pdf    Simple DIM alone
#   sim_clustered.pdf Simple DIM + Cluster-randomized
#   sim_blocked.pdf   Simple DIM + Cluster-randomized + Blocked
#   sim_adjusted.pdf  All four (adds Lin-adjusted)
# Each panel overlays on the same x- and y-axis for visual progression.

set.seed(42)
suppressPackageStartupMessages(library(sandwich))

out_dir <- path.expand("~/teaching/stat256-2026/for2026/slides/02POM/images")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- DGP (matches PS200C sim_experiments.R) -------------------------------
N_precincts <- 100
hh_per_prec <- 10
ppl_per_hh  <- 4
N           <- N_precincts * hh_per_prec * ppl_per_hh   # 4000

sigma_prec <- 1.20
sigma_hh   <- 0.70
beta_pt    <- 3.5
beta_age   <- 0.4
beta_part  <- 0.1
intercept  <- -0.3

tau_prob <- 0.05
R        <- 2000

prec_id <- rep(1:N_precincts, each = hh_per_prec * ppl_per_hh)
hh_id   <- rep(1:(N_precincts * hh_per_prec), each = ppl_per_hh)

past_turnout <- rnorm(N)
age          <- rnorm(N)
partisan     <- rnorm(N)
u_prec       <- rnorm(N_precincts, sd = sigma_prec)[prec_id]
u_hh         <- rnorm(max(hh_id),  sd = sigma_hh)[hh_id]

linpred0 <- intercept + beta_pt*past_turnout + beta_age*age +
            beta_part*partisan + u_prec + u_hh
p0 <- plogis(linpred0)
p1 <- pmin(pmax(p0 + tau_prob, 0), 1)

draw_Y <- function(T_vec) rbinom(N, 1, ifelse(T_vec == 1, p1, p0))

# Block id: quartiles of past_turnout
block_id <- as.integer(cut(past_turnout,
                           breaks = quantile(past_turnout, 0:4/4),
                           include.lowest = TRUE))

assign_simple  <- function() rbinom(N, 1, 0.5)
assign_cluster <- function() rbinom(N_precincts, 1, 0.5)[prec_id]
assign_blocked <- function() {
  T <- integer(N)
  for (b in unique(block_id)) {
    idx <- which(block_id == b)
    T[idx] <- sample(rep(c(0,1), length.out = length(idx)))
  }
  T
}

# ---- Estimators -----------------------------------------------------------
est_dim <- function(Y, T) unname(coef(lm(Y ~ T))["T"])

est_lin <- function(Y, T) {
  Xc <- cbind(past_turnout, age, partisan)
  Xc <- sweep(Xc, 2, colMeans(Xc), "-")
  colnames(Xc) <- c("pt_c", "age_c", "part_c")
  df <- data.frame(Y = Y, T = T, Xc)
  fit <- lm(Y ~ T + pt_c + age_c + part_c +
                T:pt_c + T:age_c + T:part_c, data = df)
  unname(coef(fit)["T"])
}

est_blocked <- function(Y, T) unname(coef(lm(Y ~ T + factor(block_id)))["T"])

# ---- Simulation -----------------------------------------------------------
tau_simple    <- numeric(R)
tau_clustered <- numeric(R)
tau_blocked   <- numeric(R)
tau_adjusted  <- numeric(R)

cat("Running", R, "replications...\n")
for (r in seq_len(R)) {
  T1 <- assign_simple();  Y1 <- draw_Y(T1)
  tau_simple[r]   <- est_dim(Y1, T1)
  tau_adjusted[r] <- est_lin(Y1, T1)

  T2 <- assign_blocked(); Y2 <- draw_Y(T2)
  tau_blocked[r]  <- est_blocked(Y2, T2)

  T3 <- assign_cluster(); Y3 <- draw_Y(T3)
  tau_clustered[r] <- est_dim(Y3, T3)
}
cat("Done.\n\n")

cat(sprintf("Simple DIM:         SD(actual) = %.4f\n", sd(tau_simple)))
cat(sprintf("Cluster-randomized: SD(actual) = %.4f\n", sd(tau_clustered)))
cat(sprintf("Blocked (block FE): SD(actual) = %.4f\n", sd(tau_blocked)))
cat(sprintf("Lin-adjusted:       SD(actual) = %.4f\n", sd(tau_adjusted)))

# ---- Plot setup -----------------------------------------------------------
xlim_all <- range(c(tau_simple, tau_clustered, tau_blocked, tau_adjusted))
brks     <- seq(xlim_all[1], xlim_all[2], length.out = 60)
ylim_max <- max(sapply(list(tau_simple, tau_clustered, tau_blocked, tau_adjusted),
                       function(x) max(density(x)$y))) * 1.05

col_simple_bar <- adjustcolor("grey70", alpha.f = 0.7)
col_simple     <- "grey40"
col_clustered  <- "firebrick"
col_blocked    <- "seagreen"
col_adjusted   <- "darkorange"

# Helper: build a panel with Simple DIM as bars+line and extra curves on top
draw_panel <- function(layers, fname, main) {
  pdf(file.path(out_dir, fname), width = 5.5, height = 3.4)
  par(mar = c(4, 4, 3, 1))
  hist(tau_simple, breaks = brks, col = col_simple_bar, border = "white",
       xlim = xlim_all, ylim = c(0, ylim_max), freq = FALSE,
       main = main, xlab = expression(hat(tau)))
  abline(v = tau_prob, col = "black", lwd = 2, lty = 2)
  lines(density(tau_simple), col = col_simple, lwd = 2)
  lab <- "Simple DIM"; col <- col_simple
  for (lyr in layers) {
    lines(density(lyr$x), col = lyr$col, lwd = 2.5)
    lab <- c(lab, lyr$label); col <- c(col, lyr$col)
  }
  legend("topright", legend = lab, col = col, lwd = 2.5, bty = "n", cex = 0.85)
  dev.off()
}

# ---- Panels ---------------------------------------------------------------
# 1. Simple DIM alone
draw_panel(list(), "sim_simple.pdf",
           "Simple randomization, diff-in-means")

# 2. + Cluster-randomized
draw_panel(list(list(x = tau_clustered, col = col_clustered,
                     label = "Cluster-randomized")),
           "sim_clustered.pdf",
           "Cluster-randomize at precinct level")

# 3. + Blocked
draw_panel(list(list(x = tau_clustered, col = col_clustered,
                     label = "Cluster-randomized"),
                list(x = tau_blocked, col = col_blocked,
                     label = "Blocked (block FEs)")),
           "sim_blocked.pdf",
           "Add blocking on past-turnout quartiles")

# 4. + Lin-adjusted (all four)
draw_panel(list(list(x = tau_clustered, col = col_clustered,
                     label = "Cluster-randomized"),
                list(x = tau_blocked, col = col_blocked,
                     label = "Blocked (block FEs)"),
                list(x = tau_adjusted, col = col_adjusted,
                     label = "Lin-adjusted")),
           "sim_adjusted.pdf",
           "All four: + Lin-adjusted")

cat("\nWrote figures to:", out_dir, "\n")
