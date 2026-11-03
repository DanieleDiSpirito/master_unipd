# ==============================================================================
# HOMEWORK 1 - PROBLEM 1.3 (Detailed Beginner-Friendly Version)
# ==============================================================================
# Objective:
# The de Moivre-Laplace theorem (the central limit theorem for binomial counts)
# states that as the number of trials n increases, the standardized binomial
# variable W converges in distribution to a standard normal distribution:
#
#   W = (S - n * theta) / sqrt(n * theta * (1 - theta)) ---> N(0, 1)
#
# If W becomes approximately N(0, 1), then the probability that W falls in
# ANY fixed interval [a, b], say [-0.5, -0.1], must converge to:
#
#   P(-0.5 <= W <= -0.1) ---> Phi(-0.1) - Phi(-0.5)
#
# In this exercise:
# 1. We compute the theoretical normal target value using pnorm().
# 2. For n in {30, 50, 100, 150, 300, 500} and theta in {0.1, 0.5, 0.9},
#    we simulate N = 100,000 values of W and compute the empirical fraction p_hat.
# 3. We also calculate the exact binomial probability to prove that the jumps
#    and oscillations are a real mathematical property of discrete distributions,
#    not just random noise!
# ==============================================================================

set.seed(47)

ns <- c(30, 50, 100, 150, 300, 500)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5

# The interval endpoints for W
a <- -0.5
b <- -0.1

# ------------------------------------------------------------------------------
# STEP 1: Compute the asymptotic target value from the Normal distribution
# ------------------------------------------------------------------------------
# pnorm(x) calculates the cumulative probability Phi(x) = P(Z <= x).
# P(a <= Z <= b) = P(Z <= b) - P(Z <= a) = Phi(b) - Phi(a).
target_norm <- pnorm(b) - pnorm(a)
cat("Target value from de Moivre-Laplace theorem: Phi(-0.1) - Phi(-0.5) =",
    round(target_norm, 5), "\n\n")

# ------------------------------------------------------------------------------
# STEP 2: Pre-allocate matrices to store results
# ------------------------------------------------------------------------------
p_hat   <- matrix(NA, nrow = length(ns), ncol = length(thetas),
                  dimnames = list(paste0("n=", ns), paste0("theta=", thetas)))
p_exact <- p_hat

cat("Running simulations for all (n, theta) combinations...\n")

for (j in seq_along(thetas)) {
  theta <- thetas[j]
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    
    # 1. Simulation approach (Monte Carlo):
    # Standard deviation of the binomial count S
    sigma <- sqrt(n * theta * (1 - theta))
    
    # Generate N binomial counts
    S <- rbinom(N, size = n, prob = theta)
    
    # Standardize to obtain W
    W <- (S - n * theta) / sigma
    
    # Fraction of Ws landing between -0.5 and -0.1
    p_hat[i, j] <- mean(W >= a & W <= b)
    
    # 2. Exact theoretical probability:
    # Notice that W is between a and b if and only if S satisfies:
    #   n * theta + a * sigma <= S <= n * theta + b * sigma
    # Since S can only take integer values 0, 1, ..., n, we can simply
    # sum dbinom(s, n, theta) for all integer values of s that satisfy this!
    s_all <- 0:n
    w_all <- (s_all - n * theta) / sigma
    qualifying_s <- s_all[w_all >= a & w_all <= b]
    p_exact[i, j] <- sum(dbinom(qualifying_s, size = n, prob = theta))
  }
}

cat("Simulations completed!\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Display results in a clear comparison table
# ------------------------------------------------------------------------------
cat("========================================================================\n")
cat("TABLE: CONVERGENCE OF P(-0.5 <= W <= -0.1)\n")
cat("Target Normal Limit = ", round(target_norm, 4), "\n")
cat("========================================================================\n")

for (j in seq_along(thetas)) {
  cat(sprintf("\n--- Results for theta = %0.1f ---\n", thetas[j]))
  df <- data.frame(
    n             = ns,
    Simulated_Hat = round(p_hat[, j], 4),
    Exact_Binom   = round(p_exact[, j], 4),
    Normal_Target = round(target_norm, 4),
    Absolute_Error = round(abs(p_hat[, j] - target_norm), 4)
  )
  print(df, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# STEP 4: Explaining the curious case of n = 30, theta = 0.1
# ------------------------------------------------------------------------------
# Why did we get 0.0000 for n = 30 and theta = 0.1? Let's inspect the math!
n_demo <- 30
th_demo <- 0.1
mean_demo <- n_demo * th_demo # 3
sd_demo <- sqrt(n_demo * th_demo * (1 - th_demo)) # sqrt(2.7) = 1.643168

s_lower <- mean_demo + a * sd_demo
s_upper <- mean_demo + b * sd_demo

cat("\n========================================================================\n")
cat("DEEP DIVE: WHY IS P(-0.5 <= W <= -0.1) = 0 FOR n = 30, theta = 0.1?\n")
cat("========================================================================\n")
cat(sprintf("For W to be in [%0.1f, %0.1f], S must be between:\n", a, b))
cat(sprintf("  Lower bound for S = 30*0.1 + (-0.5)*%0.4f = %0.4f\n", sd_demo, s_lower))
cat(sprintf("  Upper bound for S = 30*0.1 + (-0.1)*%0.4f = %0.4f\n", sd_demo, s_upper))
cat(sprintf("The required range for S is: [%0.4f, %0.4f]\n", s_lower, s_upper))
cat("Because S is a count, it can ONLY be an integer (..., 1, 2, 3, 4, ...).\n")
cat("There are NO integers between 2.1784 and 2.8357!\n")
cat("Therefore, it is PHYSICALLY IMPOSSIBLE for W to be in [-0.5, -0.1] at n = 30!\n")
cat("This is why the simulated and exact probabilities are both EXACTLY ZERO.\n\n")

# ------------------------------------------------------------------------------
# STEP 5: Plotting the results (2x2 layout, log scale on x-axis)
# ------------------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

for (j in seq_along(thetas)) {
  # We use log="x" because n ranges from 30 to 500 (spacing points nicely)
  plot(ns, p_hat[, j], type = "b", pch = 16, col = "blue", lwd = 2,
       ylim = c(0, 0.22), log = "x", xaxt = "n",
       xlab = "Sample size n (logarithmic scale)",
       ylab = expression(hat(p)[list(n, theta)]),
       main = bquote("Convergence for" ~ theta == .(thetas[j])))
  
  # Custom x-axis showing actual n values
  axis(1, at = ns, labels = ns)
  
  # Overlay exact binomial points as red crosses
  points(ns, p_exact[, j], pch = 4, col = "red", cex = 1.3, lwd = 2)
  
  # Horizontal dashed line for the asymptotic normal limit
  abline(h = target_norm, lty = 2, col = "black", lwd = 1.5)
}

plot.new()
legend("center",
       legend = c("Simulated frequency (p_hat)", "Exact binomial probability", "Normal limit (Phi(-0.1) - Phi(-0.5))"),
       col = c("blue", "red", "black"),
       pch = c(16, 4, NA), lty = c(1, NA, 2), lwd = 2,
       bty = "n", cex = 1.05, title = "Legend")
