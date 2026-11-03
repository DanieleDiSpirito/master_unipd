# HW1 - Problem 1.1
# coverage of exact, Wald and score 0.90 CIs for the binomial parameter

set.seed(47)
ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5
alpha <- 0.10
z <- qnorm(1 - alpha/2)

# quick check: the beta formula gives the same interval as binom.test
# and the Wilson formula the same as prop.test(correct = FALSE)
s <- 3; n <- 10; p <- s/n
binom.test(s, n, conf.level = 0.9)$conf.int
c(qbeta(alpha/2, s, n - s + 1), qbeta(1 - alpha/2, s + 1, n - s))
prop.test(s, n, conf.level = 0.9, correct = FALSE)$conf.int
(p + z^2/(2*n) + c(-1, 1) * z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)

# matrices for the results (rows = n, columns = theta)
cov_exact <- matrix(NA, length(ns), length(thetas), dimnames = list(ns, thetas))
cov_wald  <- cov_exact
cov_score <- cov_exact

for (i in 1:length(ns)) {
  for (j in 1:length(thetas)) {
    n <- ns[i]
    theta <- thetas[j]
    S <- rbinom(N, n, theta)
    p <- S/n

    # exact (Clopper-Pearson)
    lo <- qbeta(alpha/2, S, n - S + 1)
    up <- qbeta(1 - alpha/2, S + 1, n - S)
    cov_exact[i, j] <- mean(lo <= theta & theta <= up)

    # Wald (slides)
    lo <- p - z * sqrt(p*(1-p)/n)
    up <- p + z * sqrt(p*(1-p)/n)
    cov_wald[i, j] <- mean(lo <= theta & theta <= up)

    # score / Wilson (prop.test with correct = FALSE)
    lo <- (p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
    up <- (p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
    cov_score[i, j] <- mean(lo <= theta & theta <= up)
  }
}

round(cov_exact, 4)
round(cov_wald, 4)
round(cov_score, 4)

# plot: one panel for each theta, the 4th panel is the legend
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in 1:length(thetas)) {
  matplot(ns, cbind(cov_exact[, j], cov_wald[, j], cov_score[, j]),
          type = "b", pch = 15:17, lty = 1, lwd = 2,
          col = c("black", "red", "blue"), ylim = c(0.6, 1),
          xlab = "n", ylab = "coverage", main = paste("theta =", thetas[j]))
  abline(h = 0.9, lty = 2)
}
plot.new()
legend("center", c("exact", "Wald", "score"), col = c("black", "red", "blue"),
       pch = 15:17, lwd = 2, bty = "n")
