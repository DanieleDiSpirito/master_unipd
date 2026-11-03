# Problem 1.1
# Exact, Wald, and Score confidence intervals using binom.test and prop.test directly

set.seed(47)

ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5
alpha <- 0.10
z <- qnorm(1 - alpha/2)

results <- data.frame()

for (n in ns) {
  # 1. Exact interval using binom.test directly
  exact_ci <- sapply(0:n, function(s) {
    binom.test(s, n, conf.level = 0.90)$conf.int
  })
  
  # 2. Wald interval using the Z formula from the slides
  wald_ci <- sapply(0:n, function(s) {
    p <- s / n
    p + c(-1, 1) * z * sqrt(p * (1 - p) / n)
  })
  
  # 3. Approximate interval using prop.test with correct = FALSE directly
  score_ci <- sapply(0:n, function(s) {
    prop.test(s, n, conf.level = 0.90, correct = FALSE)$conf.int
  })
  
  for (theta in thetas) {
    # Generate N binomial random numbers
    S <- rbinom(N, size = n, prob = theta)
    
    # Check the fraction of times theta is inside each interval
    cov_exact <- mean(exact_ci[1, S + 1] <= theta & theta <= exact_ci[2, S + 1])
    cov_wald  <- mean(wald_ci[1, S + 1] <= theta & theta <= wald_ci[2, S + 1])
    cov_score <- mean(score_ci[1, S + 1] <= theta & theta <= score_ci[2, S + 1])
    
    results <- rbind(results, data.frame(
      n = n, theta = theta,
      Exact = cov_exact, Wald = cov_wald, Score = cov_score
    ))
  }
}

print(round(results, 4))

# Plot results
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (th in thetas) {
  sub_res <- results[results$theta == th, ]
  matplot(sub_res$n, sub_res[, c("Exact", "Wald", "Score")],
          type = "b", pch = 15:17, lty = 1, lwd = 2,
          col = c("black", "red", "blue"), ylim = c(0.6, 1.0),
          xlab = "n", ylab = "Coverage", main = paste("theta =", th))
  abline(h = 0.90, lty = 2, col = "gray50")
}
plot.new()
legend("center", legend = c("Exact (binom.test)", "Wald (slides)", "Score (prop.test)"),
       col = c("black", "red", "blue"), pch = 15:17, lwd = 2, bty = "n")
