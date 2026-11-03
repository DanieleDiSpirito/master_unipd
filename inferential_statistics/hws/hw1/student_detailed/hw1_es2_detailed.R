# ==============================================================================
# HOMEWORK 1 - PROBLEM 1.2 (Detailed Beginner-Friendly Version)
# ==============================================================================
# Objective:
# Investigate the effect of the CONTINUITY CORRECTION (CC) on binomial
# confidence intervals, following the classic paper by Newcombe (1998).
#
# Why do we need a continuity correction?
# Because the binomial distribution is DISCRETE (counts S can only be 0, 1, 2...)
# whereas the normal distribution we use as an approximation is CONTINUOUS.
# A single discrete count s is represented in the continuous world by a bar
# running from (s - 0.5) to (s + 0.5).
# Dividing by n, the correction shifts the sample proportion by 0.5 / n = 1 / (2n).
#
# We compare 6 methods:
#   1. Wald (uncorrected)
#   2. Wald + CC: pushes bounds outward by 1/(2n)
#   3. Wilson Score (uncorrected): prop.test(correct = FALSE)
#   4. Wilson Score + CC: prop.test(correct = TRUE)
#   5. Exact Mid-p interval: the "uncorrected" version of Clopper-Pearson
#   6. Exact Clopper-Pearson: binom.test(), which already includes full endpoint
# ==============================================================================

set.seed(47)

ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5
conf_level <- 0.90
alpha <- 1 - conf_level
z <- qnorm(1 - alpha / 2) # 1.644854

# ------------------------------------------------------------------------------
# HELPER FUNCTION 1: Wilson Score Formula
# ------------------------------------------------------------------------------
# This function calculates the Wilson score interval for any proportion vector p
# and sample size n.
wilson_formula <- function(p, n) {
  denom  <- 1 + (z^2) / n
  center <- (p + (z^2) / (2 * n)) / denom
  margin <- (z * sqrt((p * (1 - p) / n) + (z^2) / (4 * n^2))) / denom
  cbind(center - margin, center + margin)
}

# ------------------------------------------------------------------------------
# HELPER FUNCTION 2: Mid-p Exact Interval
# ------------------------------------------------------------------------------
# In Clopper-Pearson, the tail probability includes the full probability P(S = s).
# In Mid-p, we count only HALF of P(S = s).
# Mathematically:
#   Lower bound solves: P_theta(S > s) + 0.5 * P_theta(S = s) = alpha / 2
#   Upper bound solves: P_theta(S < s) + 0.5 * P_theta(S = s) = alpha / 2
# Since there is no closed-form algebraic formula, we solve it numerically
# using R's built-in 1-dimensional solver: uniroot().
compute_midp_for_s <- function(s, n) {
  lo <- 0
  up <- 1
  
  if (s > 0) {
    # 1 - pbinom(s, n, th) is P(S > s). We subtract alpha/2 to find where it equals 0.
    f_lower <- function(th) {
      (1 - pbinom(s, n, th)) + 0.5 * dbinom(s, n, th) - (alpha / 2)
    }
    # uniroot searches for the value of th in [0, 1] where f_lower(th) == 0
    lo <- uniroot(f_lower, interval = c(0, 1), tol = 1e-8)$root
  }
  
  if (s < n) {
    # pbinom(s - 1, n, th) is P(S < s).
    f_upper <- function(th) {
      pbinom(s - 1, n, th) + 0.5 * dbinom(s, n, th) - (alpha / 2)
    }
    up <- uniroot(f_upper, interval = c(0, 1), tol = 1e-8)$root
  }
  
  return(c(lo, up))
}

# ------------------------------------------------------------------------------
# Arrays to store coverage and interval width for each method
# ------------------------------------------------------------------------------
method_names <- c("Wald", "Wald+CC", "Score", "Score+CC", "Mid-p", "Exact (C-P)")

# 3D arrays: dimensions are [sample size n, theta, method]
cov_array   <- array(NA, dim = c(length(ns), length(thetas), length(method_names)),
                     dimnames = list(paste0("n=", ns), paste0("th=", thetas), method_names))
width_array <- cov_array

cat("Starting simulation for Problem 1.2 (Coverage and Width)...\n")

for (i in seq_along(ns)) {
  n <- ns[i]
  
  # Optimization tip for newbies:
  # Since s can only take integer values 0, 1, ..., n, the mid-p interval can only
  # have (n + 1) possible values! Instead of running uniroot 100,000 times,
  # we precompute it once for each s in 0:n and look it up in a table!
  midp_lookup <- t(sapply(0:n, function(s) compute_midp_for_s(s, n)))
  
  for (j in seq_along(thetas)) {
    theta <- thetas[j]
    
    # Simulate N experiments
    S <- rbinom(N, size = n, prob = theta)
    p <- S / n
    se <- sqrt(p * (1 - p) / n)
    
    # Store lower and upper limits in named lists
    lo <- list()
    up <- list()
    
    # 1. Wald (uncorrected)
    lo[["Wald"]] <- p - z * se
    up[["Wald"]] <- p + z * se
    
    # 2. Wald + CC: widen by adding 1/(2n) to the margin of error
    # pmax(0, ...) and pmin(1, ...) ensure the interval doesn't go below 0 or above 1
    lo[["Wald+CC"]] <- pmax(0, p - z * se - 1 / (2 * n))
    up[["Wald+CC"]] <- pmin(1, p + z * se + 1 / (2 * n))
    
    # 3. Score (uncorrected Wilson)
    w_uncorr <- wilson_formula(p, n)
    lo[["Score"]] <- w_uncorr[, 1]
    up[["Score"]] <- w_uncorr[, 2]
    
    # 4. Score + CC (as implemented by prop.test(correct = TRUE))
    # Shifts p down by 0.5/n for lower limit, and up by 0.5/n for upper limit.
    # When S = n/2 (exact middle), no shift is applied.
    cc_shift <- ifelse(S == n / 2, 0, 0.5 / n)
    lo[["Score+CC"]] <- ifelse(p - cc_shift <= 0, 0, wilson_formula(p - cc_shift, n)[, 1])
    up[["Score+CC"]] <- ifelse(p + cc_shift >= 1, 1, wilson_formula(p + cc_shift, n)[, 2])
    
    # 5. Mid-p exact
    # We look up the precomputed table: S = 0 corresponds to row 1, so S + 1
    lo[["Mid-p"]] <- midp_lookup[S + 1, 1]
    up[["Mid-p"]] <- midp_lookup[S + 1, 2]
    
    # 6. Clopper-Pearson exact (via Beta distribution quantiles)
    lo[["Exact (C-P)"]] <- qbeta(alpha / 2, S, n - S + 1)
    up[["Exact (C-P)"]] <- qbeta(1 - alpha / 2, S + 1, n - S)
    
    # Evaluate empirical coverage and average width for each method
    for (m in method_names) {
      cov_array[i, j, m]   <- mean(lo[[m]] <= theta & theta <= up[[m]])
      width_array[i, j, m] <- mean(up[[m]] - lo[[m]])
    }
  }
}

cat("Simulation completed successfully!\n\n")

# ------------------------------------------------------------------------------
# STEP 4: Reporting results in tables
# ------------------------------------------------------------------------------
cat("========================================================================\n")
cat("PART 1: EMPIRICAL COVERAGE (Target = 0.90)\n")
cat("========================================================================\n")
for (j in seq_along(thetas)) {
  cat(sprintf("\n--- Coverage for True theta = %0.1f ---\n", thetas[j]))
  print(round(cov_array[, j, ], 4))
}

cat("\n========================================================================\n")
cat("PART 2: AVERAGE INTERVAL WIDTH (Smaller is narrower/more informative)\n")
cat("========================================================================\n")
for (j in seq_along(thetas)) {
  cat(sprintf("\n--- Average Width for True theta = %0.1f ---\n", thetas[j]))
  print(round(width_array[, j, ], 3))
}

# ------------------------------------------------------------------------------
# STEP 5: Visualizing Coverage and Width (2x2 layout)
# ------------------------------------------------------------------------------
# Solid lines = without CC; Dashed lines = with CC
cols <- c("red", "red", "blue", "blue", "darkgreen", "darkgreen")
ltys <- c(1, 2, 1, 2, 1, 2)
pchs <- c(16, 1, 17, 2, 15, 0)

# FIGURE 1: COVERAGE
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
for (j in seq_along(thetas)) {
  matplot(ns, cov_array[, j, ], type = "b", pch = pchs, col = cols, lty = ltys,
          lwd = 2, ylim = c(0.6, 1.0),
          xlab = "Sample size (n)", ylab = "Empirical Coverage Fraction",
          main = bquote("Coverage:" ~ theta == .(thetas[j])))
  abline(h = conf_level, lty = 3, col = "gray40", lwd = 1.5)
}
plot.new()
legend("center", legend = method_names, col = cols, lty = ltys, pch = pchs,
       lwd = 2, bty = "n", cex = 1.05, title = "Method")

# FIGURE 2: WIDTH
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
for (j in seq_along(thetas)) {
  matplot(ns, width_array[, j, ], type = "b", pch = pchs, col = cols, lty = ltys,
          lwd = 2, ylim = c(0, 0.6),
          xlab = "Sample size (n)", ylab = "Average Interval Width",
          main = bquote("Width:" ~ theta == .(thetas[j])))
}
plot.new()
legend("center", legend = method_names, col = cols, lty = ltys, pch = pchs,
       lwd = 2, bty = "n", cex = 1.05, title = "Method")
