# Problem 1.3
# Checking de Moivre-Laplace theorem: W -> N(0, 1)

set.seed(47)

ns <- c(30, 50, 100, 150, 300, 500)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5

# Target value from the normal distribution
target <- pnorm(-0.1) - pnorm(-0.5)
cat("Normal target:", round(target, 4), "\n")

res <- data.frame()

for (theta in thetas) {
  for (n in ns) {
    sigma <- sqrt(n * theta * (1 - theta))
    
    # Simulation
    S <- rbinom(N, size = n, prob = theta)
    W <- (S - n * theta) / sigma
    p_hat <- mean(W >= -0.5 & W <= -0.1)
    
    # Exact binomial probability
    s <- 0:n
    w <- (s - n * theta) / sigma
    p_exact <- sum(dbinom(s[w >= -0.5 & w <= -0.1], n, theta))
    
    res <- rbind(res, data.frame(
      theta = theta, n = n,
      Simulated = p_hat, Exact = p_exact, Target = target
    ))
  }
}

print(round(res, 4))

# Plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (th in thetas) {
  sub_r <- res[res$theta == th, ]
  plot(sub_r$n, sub_r$Simulated, type = "b", pch = 16, col = "blue", lwd = 2,
       ylim = c(0, 0.22), log = "x", xlab = "n (log scale)", ylab = "p_hat",
       main = paste("theta =", th))
  points(sub_r$n, sub_r$Exact, pch = 4, col = "red", cex = 1.3, lwd = 2)
  abline(h = target, lty = 2)
}
plot.new()
legend("center", c("Simulated", "Exact", "Target"),
       col = c("blue", "red", "black"), pch = c(16, 4, NA), lty = c(1, NA, 2),
       lwd = 2, bty = "n")
