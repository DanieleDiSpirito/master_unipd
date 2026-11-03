# ==============================================================================
# HOMEWORK 1 - PROBLEM 1.1 (Detailed Beginner-Friendly Version)
# ==============================================================================
# Objective:
# Compare three different methods for constructing a 90% confidence interval
# for a binomial success probability theta:
#   1. Exact (Clopper-Pearson) interval: from binom.test()
#   2. Wald interval (approximate Z formula from course slides)
#   3. Wilson Score interval: from prop.test(..., correct = FALSE)
#
# We test these methods across different sample sizes n in {10, 30, 50, 100}
# and true probabilities theta in {0.1, 0.5, 0.9}.
# For each setting, we generate N = 100,000 datasets and measure what fraction
# of the computed intervals actually "trap" (contain) the true theta.
# ==============================================================================

# ------------------------------------------------------------------------------
# STEP 1: Setting up the simulation parameters
# ------------------------------------------------------------------------------
# Setting a random seed ensures that anyone running this script gets the EXACT
# same random numbers and results.
set.seed(47)

# The sample sizes n we are asked to investigate
ns <- c(10, 30, 50, 100)

# The true values of the success parameter theta
thetas <- c(0.1, 0.5, 0.9)

# Number of simulated experiments for each combination of (n, theta).
# 100,000 is large enough so that random simulation noise is less than 0.1%.
N <- 1e5

# Confidence level is 90%, which means the error rate alpha is 10% (0.10).
# For a two-sided interval, alpha is split equally between the two tails: 5% each.
conf_level <- 0.90
alpha <- 1 - conf_level

# z is the critical value from the standard normal distribution N(0, 1).
# For 90% confidence, we need the 95th percentile (1 - alpha/2 = 0.95).
# qnorm(0.95) is approximately 1.644854.
z <- qnorm(1 - alpha / 2)

cat("Critical value z for 90% confidence interval:", round(z, 4), "\n\n")

# ------------------------------------------------------------------------------
# STEP 2: Pre-allocating matrices to store our empirical results
# ------------------------------------------------------------------------------
# Pre-allocating matrices with known dimensions is much faster and cleaner
# than growing vectors inside loops.
# Rows will represent different sample sizes (n), columns represent thetas.
cov_exact <- matrix(NA, nrow = length(ns), ncol = length(thetas),
                    dimnames = list(paste0("n=", ns), paste0("theta=", thetas)))
cov_wald  <- cov_exact
cov_score <- cov_exact

# ------------------------------------------------------------------------------
# STEP 3: Running the simulation loop
# ------------------------------------------------------------------------------
cat("Starting simulation over all (n, theta) combinations...\n")

for (i in seq_along(ns)) {
  n <- ns[i]
  
  for (j in seq_along(thetas)) {
    theta <- thetas[j]
    
    # rbinom(N, size = n, prob = theta) simulates N independent experiments.
    # In each experiment, we flip an unfair coin with success probability theta
    # n times, and S stores the number of successes (an integer between 0 and n).
    S <- rbinom(N, size = n, prob = theta)
    
    # p is our sample proportion of successes: hat(theta) = S / n.
    p <- S / n
    
    # --------------------------------------------------------------------------
    # Method 1: Exact Clopper-Pearson Interval
    # --------------------------------------------------------------------------
    # binom.test(s, n)$conf.int finds the limits by inverting the binomial tails.
    # In probability theory, there is a known exact mathematical equivalence:
    # the lower limit corresponds to the alpha/2 quantile of a Beta(s, n - s + 1),
    # and the upper limit is the (1 - alpha/2) quantile of a Beta(s + 1, n - s).
    # Using qbeta() on the vector S computes all 100,000 intervals in milliseconds,
    # giving the EXACT same numbers as calling binom.test() 100,000 times!
    lo_exact <- qbeta(alpha / 2, S, n - S + 1)
    up_exact <- qbeta(1 - alpha / 2, S + 1, n - S)
    
    # We check if the true theta falls inside [lo, up].
    # In R, (lo <= theta & theta <= up) returns a vector of TRUE / FALSE.
    # Taking the mean() of a logical vector computes the proportion of TRUEs.
    cov_exact[i, j] <- mean(lo_exact <= theta & theta <= up_exact)
    
    # --------------------------------------------------------------------------
    # Method 2: Wald Approximate Interval (Z-formula from slides)
    # --------------------------------------------------------------------------
    # Formula: p +- z * sqrt(p * (1 - p) / n)
    # Notice that it estimates the standard error using the sample proportion p!
    # When S = 0, p = 0, so the standard error is sqrt(0) = 0.
    # In that case, the interval collapses to a single point [0, 0]!
    se_wald <- sqrt(p * (1 - p) / n)
    lo_wald <- p - z * se_wald
    up_wald <- p + z * se_wald
    
    cov_wald[i, j] <- mean(lo_wald <= theta & theta <= up_wald)
    
    # --------------------------------------------------------------------------
    # Method 3: Wilson Score Interval (prop.test with correct = FALSE)
    # --------------------------------------------------------------------------
    # Instead of plugging in p into the standard error, Wilson (1927) tests
    # H0: theta = theta0 with the TRUE variance theta0*(1 - theta0)/n and
    # inverts the resulting quadratic inequality.
    # This leads to a formula with an adjusted center and variance:
    denom  <- 1 + (z^2) / n
    center <- (p + (z^2) / (2 * n)) / denom
    margin <- (z * sqrt((p * (1 - p) / n) + (z^2) / (4 * n^2))) / denom
    lo_score <- center - margin
    up_score <- center + margin
    
    cov_score[i, j] <- mean(lo_score <= theta & theta <= up_score)
  }
}

cat("Simulation completed successfully!\n\n")

# ------------------------------------------------------------------------------
# STEP 4: Displaying the results in clear tables
# ------------------------------------------------------------------------------
cat("========================================================================\n")
cat("TABLE: EMPIRICAL COVERAGE (Nominal Level = 0.90)\n")
cat("========================================================================\n")

for (j in seq_along(thetas)) {
  th <- thetas[j]
  cat(sprintf("\n--- Results for True theta = %0.1f ---\n", th))
  tab <- data.frame(
    n = ns,
    Exact_Clopper_Pearson = round(cov_exact[, j], 4),
    Wald_Z_Formula        = round(cov_wald[, j], 4),
    Wilson_Score          = round(cov_score[, j], 4)
  )
  print(tab, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# STEP 5: Visualizing the results (2x2 layout so no plot is cut off)
# ------------------------------------------------------------------------------
# We divide the plotting window into 2 rows and 2 columns:
# - Panel 1 (top-left): theta = 0.1
# - Panel 2 (top-right): theta = 0.5
# - Panel 3 (bottom-left): theta = 0.9
# - Panel 4 (bottom-right): clean legend placed in the empty 4th box!
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

cols <- c("black", "red", "blue")
pchs <- c(15, 16, 17) # square, circle, triangle

for (j in seq_along(thetas)) {
  # matplot plots all 3 methods on the same graph against n
  matplot(ns, cbind(cov_exact[, j], cov_wald[, j], cov_score[, j]),
          type = "b", pch = pchs, lty = 1, lwd = 2, col = cols,
          ylim = c(0.6, 1.0),
          xlab = "Sample size (n)",
          ylab = "Empirical Coverage Fraction",
          main = bquote("True" ~ theta == .(thetas[j])))
  
  # Add a horizontal dashed line at the nominal 90% target
  abline(h = conf_level, lty = 2, col = "gray40", lwd = 1.5)
}

# 4th panel: blank canvas reserved exclusively for the legend
plot.new()
legend("center",
       legend = c("Exact (Clopper-Pearson)", "Wald (Z-formula)", "Wilson Score (prop.test)"),
       col = cols, pch = pchs, lty = 1, lwd = 2,
       bty = "n", cex = 1.1, title = "Confidence Interval Method")
