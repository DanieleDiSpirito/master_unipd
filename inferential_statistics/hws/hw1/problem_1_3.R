# Problem 1.3 -- de Moivre-Laplace theorem by simulation.
#
# S ~ Bin(n, theta) and W = (S - n theta) / sqrt(n theta (1 - theta)).
# By de Moivre-Laplace, W -> N(0, 1) in distribution, hence
#   p_{n,theta} = P(-0.5 <= W <= -0.1) -> Phi(-0.1) - Phi(-0.5).
# For each theta and n we approximate p_{n,theta} by the relative frequency
# hat p_{n,theta} over N = 1e5 simulated values of W.
# For reference we also compute the exact value of p_{n,theta} from the
# binomial distribution.

set.seed(47)
ns     <- c(30, 50, 100, 150, 300, 500)
thetas <- c(0.1, 0.5, 0.9)
N      <- 1e5
a <- -0.5
b <- -0.1

target <- pnorm(b) - pnorm(a)    # limit value given by the theorem

res <- expand.grid(n = ns, theta = thetas)
res$p_hat   <- NA_real_
res$p_exact <- NA_real_

for (i in seq_len(nrow(res))) {
  n     <- res$n[i]
  theta <- res$theta[i]
  sd    <- sqrt(n * theta * (1 - theta))
  S <- rbinom(N, size = n, prob = theta)
  W <- (S - n * theta) / sd
  res$p_hat[i] <- mean(a <= W & W <= b)
  # exact: P(n theta + a sd <= S <= n theta + b sd)
  s_vals <- 0:n
  w_vals <- (s_vals - n * theta) / sd
  res$p_exact[i] <- sum(dbinom(s_vals[a <= w_vals & w_vals <= b], n, theta))
}

cat("Target value Phi(-0.1) - Phi(-0.5) =", round(target, 4), "\n")
print(res, digits = 4)

# ---- plot ------------------------------------------------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
ylim <- range(c(0, res$p_hat, res$p_exact, target))
for (theta in thetas) {
  d <- res[res$theta == theta, ]
  plot(d$n, d$p_hat, type = "b", pch = 16, lwd = 2, col = "blue",
       ylim = ylim, xlab = "n", ylab = expression(hat(p)[list(n, theta)]),
       main = bquote(theta == .(theta)), log = "x", xaxt = "n")
  axis(1, at = ns)
  points(d$n, d$p_exact, pch = 4, col = "red", cex = 1.5, lwd = 2)
  abline(h = target, lty = 2, lwd = 2)
}
plot.new()   # 4th panel: legend
legend("center", bty = "n",
       legend = c("Simulated (relative frequency)", "Exact binomial",
                  "Normal limit"),
       col = c("blue", "red", "black"), pch = c(16, 4, NA),
       lty = c(1, NA, 2), lwd = 2)
