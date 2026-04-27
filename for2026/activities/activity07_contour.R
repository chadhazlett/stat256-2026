# Activity 7: contour plot for sensitivity activity.
# Uses sensemakr's raw-input API for the contours, then adds 1x/2x/3x X_1
# benchmark markers manually for clearer pedagogy.
#
# Setup matches the augmented regression table in Q1:
#   Estimate = 0.30, SE = 0.10, df = 600
#   Strongest observed covariate X_1: R^2_{D~X_1|...} = 4%, R^2_{Y~X_1|...,D} = 8%

library(sensemakr)

estimate <- 0.30
se       <- 0.10
dof      <- 600

# X_1 benchmark stats (k=1)
r2dz.x1  <- 0.04
r2yz.x1  <- 0.08

ks       <- c(1, 2, 3)
xs       <- ks * r2dz.x1   # 0.04, 0.08, 0.12
ys       <- ks * r2yz.x1   # 0.08, 0.16, 0.24
labs     <- paste0(ks, "x X_1")

add_markers <- function() {
  points(xs, ys, pch = 23, col = "red", bg = "red", cex = 1.2)
  text(xs, ys, labs, pos = 4, cex = 0.85, col = "red")
}

pdf("activity07_contour.pdf", width = 10, height = 4.6)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

ovb_contour_plot(estimate, se, dof,
                 sensitivity.of = "estimate",
                 main           = "Adjusted point estimate")
add_markers()

ovb_contour_plot(estimate, se, dof,
                 sensitivity.of = "t-value",
                 main           = "Adjusted t-statistic")
add_markers()

dev.off()

cat("Wrote activity07_contour.pdf\n")
