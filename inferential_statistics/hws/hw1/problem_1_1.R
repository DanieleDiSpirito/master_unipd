# Problem 1.1 -- Exact vs approximate 0.90 confidence intervals for the
# parameter theta of a binomial distribution.
#
# For every n in {10, 30, 50, 100} and theta in {0.1, 0.5, 0.9} we generate
# N = 1e5 binomial counts S ~ Bin(n, theta) and, for each of them, build
#   1. the exact (Clopper-Pearson) interval           -> binom.test
#   2. the approximate Wald (Z) interval of the slides -> S/n +- z * sqrt(p(1-p)/n)
#   3. the Wilson score interval                      -> prop.test(correct = FALSE)
# and record the fraction of intervals that trap theta (empirical coverage).
#
# Implementation note: S can only take the n + 1 values 0, 1, ..., n, so each
# interval is computed once for every possible value s and then looked up for
# the N simulated counts (row s + 1). The result is identical to calling
# binom.test / prop.test on each of the 1e5 data sets, but much faster.

set.seed(47)
ns     <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N      <- 1e5
level  <- 0.90
z      <- qnorm(1 - (1 - level) / 2)   # 1.645

# ---- intervals for every possible outcome s = 0, ..., n --------------------
# each function returns a (n + 1) x 2 matrix: [lower, upper]

ci_exact <- function(n, level) {
  t(sapply(0:n, function(s) binom.test(s, n, conf.level = level)$conf.int))
}

ci_wald <- function(n, level) {
  z <- qnorm(1 - (1 - level) / 2)
  p <- (0:n) / n
  cbind(p - z * sqrt(p * (1 - p) / n), p + z * sqrt(p * (1 - p) / n))
}

ci_score <- function(n, level) {
  t(sapply(0:n, function(s)
    prop.test(s, n, conf.level = level, correct = FALSE)$conf.int))
}

methods <- list(
  "Exact (Clopper-Pearson)"      = ci_exact,
  "Wald (Z)"                     = ci_wald,
  "Score (prop.test, no corr.)"  = ci_score
)

# ---- simulation ------------------------------------------------------------
res <- expand.grid(n = ns, theta = thetas, method = names(methods),
                   stringsAsFactors = FALSE)
res$coverage <- NA_real_

for (n in ns) {
  # intervals for all possible outcomes, for every method
  cis <- lapply(methods, function(f) f(n, level))
  for (theta in thetas) {
    S <- rbinom(N, size = n, prob = theta)   # N simulated data sets
    for (m in names(methods)) {
      lo <- cis[[m]][S + 1, 1]
      up <- cis[[m]][S + 1, 2]
      res$coverage[res$n == n & res$theta == theta & res$method == m] <-
        mean(lo <= theta & theta <= up)
    }
  }
}

# ---- table -----------------------------------------------------------------
tab <- reshape(res, idvar = c("theta", "n"), timevar = "method",
               direction = "wide")
names(tab) <- sub("coverage.", "", names(tab), fixed = TRUE)
tab <- tab[order(tab$theta, tab$n), ]
rownames(tab) <- NULL
print(tab, digits = 4)

# ---- plot ------------------------------------------------------------------
cols <- c("black", "red", "blue")
pchs <- c(15, 16, 17)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (theta in thetas) {
  plot(NA, xlim = range(ns), ylim = c(0.6, 1), xaxt = "n",
       xlab = "n", ylab = "Empirical coverage",
       main = bquote(theta == .(theta)))
  axis(1, at = ns)
  abline(h = level, lty = 2, col = "gray40")
  for (k in seq_along(methods)) {
    d <- res[res$theta == theta & res$method == names(methods)[k], ]
    lines(d$n, d$coverage, type = "b", col = cols[k], pch = pchs[k], lwd = 2)
  }
}
plot.new()   # 4th panel: legend
legend("center", legend = names(methods), col = cols, pch = pchs, lwd = 2,
       bty = "n", title = "Method")
