# HW1 - Problem 1.3
# de Moivre-Laplace: W = (S - n*theta)/sqrt(n*theta*(1-theta)) -> N(0,1)
# so P(-0.5 <= W <= -0.1) -> pnorm(-0.1) - pnorm(-0.5)

set.seed(47)
ns <- c(30, 50, 100, 150, 300, 500)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5

target <- pnorm(-0.1) - pnorm(-0.5)
target

p_hat <- matrix(NA, length(ns), length(thetas), dimnames = list(ns, thetas))
p_true <- p_hat   # exact value with the binomial, just to compare

for (j in 1:length(thetas)) {
  for (i in 1:length(ns)) {
    n <- ns[i]
    theta <- thetas[j]
    S <- rbinom(N, n, theta)
    W <- (S - n*theta) / sqrt(n*theta*(1-theta))
    p_hat[i, j] <- mean(W >= -0.5 & W <= -0.1)

    s <- 0:n
    w <- (s - n*theta) / sqrt(n*theta*(1-theta))
    p_true[i, j] <- sum(dbinom(s[w >= -0.5 & w <= -0.1], n, theta))
  }
}

round(p_hat, 4)
round(p_true, 4)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in 1:length(thetas)) {
  plot(ns, p_hat[, j], type = "b", pch = 16, col = "blue", lwd = 2,
       ylim = c(0, 0.2), log = "x", xlab = "n (log scale)", ylab = "p hat",
       main = paste("theta =", thetas[j]))
  points(ns, p_true[, j], pch = 4, col = "red", cex = 1.5)
  abline(h = target, lty = 2)
}
plot.new()
legend("center", c("simulated", "exact binomial", "normal limit"),
       col = c("blue", "red", "black"), pch = c(16, 4, NA),
       lty = c(1, NA, 2), bty = "n")

# n = 30, theta = 0.1: no integer s gives W in [-0.5, -0.1]
3 + c(-0.5, -0.1) * sqrt(30*0.1*0.9)
