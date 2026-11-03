# Problem 1.2 -- Continuity correction (CC) for binomial confidence intervals
# (Newcombe, 1998).
#
# The experiment of Problem 1.1 is repeated (same seed => same simulated data
# sets) applying the continuity correction to every method:
#
#   Wald (Z)          : p +- z sqrt(p(1-p)/n)           (Newcombe method 1)
#   Wald + CC         : p +- [z sqrt(p(1-p)/n) + 1/(2n)] (Newcombe method 2)
#   Score (Wilson)    : prop.test(correct = FALSE)       (Newcombe method 3)
#   Score + CC        : prop.test(correct = TRUE)        (Newcombe method 4)
#   Exact mid-p       : Clopper-Pearson without CC       (Newcombe method 6)
#   Exact (C-P) = CC  : binom.test                       (Newcombe method 5)
#
# The Clopper-Pearson interval is the "continuity corrected" member of the
# exact family: it counts the whole probability P(S = s) in the tails.
# Its uncorrected counterpart is the mid-p interval, which only counts
# P(S = s) / 2. So, for the exact method, "uncorrected vs corrected" means
# "mid-p vs Clopper-Pearson".
#
# As in Problem 1.1, every interval is computed once for each possible outcome
# s = 0, ..., n and then looked up for the N simulated counts.

set.seed(47)
ns     <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N      <- 1e5
level  <- 0.90
alpha  <- 1 - level

# ---- intervals for every possible outcome s = 0, ..., n --------------------

ci_wald <- function(n, level, cc = FALSE) {
  z <- qnorm(1 - (1 - level) / 2)
  p <- (0:n) / n
  half <- z * sqrt(p * (1 - p) / n) + cc / (2 * n)
  cbind(pmax(0, p - half), pmin(1, p + half))   # truncated to [0, 1]
}

ci_score <- function(n, level, cc = FALSE) {
  t(sapply(0:n, function(s)
    prop.test(s, n, conf.level = level, correct = cc)$conf.int))
}

ci_exact <- function(n, level) {               # Clopper-Pearson
  t(sapply(0:n, function(s) binom.test(s, n, conf.level = level)$conf.int))
}

ci_midp <- function(n, level) {                # mid-p exact interval
  a <- (1 - level) / 2
  t(sapply(0:n, function(s) {
    # lower: P(S > s) + P(S = s)/2 = a ; upper: P(S < s) + P(S = s)/2 = a
    lo <- if (s == 0) 0 else uniroot(function(th)
      pbinom(s, n, th, lower.tail = FALSE) + 0.5 * dbinom(s, n, th) - a,
      c(0, 1), tol = 1e-10)$root
    up <- if (s == n) 1 else uniroot(function(th)
      pbinom(s - 1, n, th) + 0.5 * dbinom(s, n, th) - a,
      c(0, 1), tol = 1e-10)$root
    c(lo, up)
  }))
}

methods <- list(
  "Wald"            = function(n, l) ci_wald(n, l, cc = FALSE),
  "Wald + CC"       = function(n, l) ci_wald(n, l, cc = TRUE),
  "Score"           = function(n, l) ci_score(n, l, cc = FALSE),
  "Score + CC"      = function(n, l) ci_score(n, l, cc = TRUE),
  "Exact mid-p"     = ci_midp,
  "Exact C-P (CC)"  = ci_exact
)

# ---- simulation (same random numbers as in Problem 1.1) -------------------
res <- expand.grid(n = ns, theta = thetas, method = names(methods),
                   stringsAsFactors = FALSE)
res$coverage <- NA_real_
res$width    <- NA_real_

for (n in ns) {
  cis <- lapply(methods, function(f) f(n, level))
  for (theta in thetas) {
    S <- rbinom(N, size = n, prob = theta)
    for (m in names(methods)) {
      lo <- cis[[m]][S + 1, 1]
      up <- cis[[m]][S + 1, 2]
      idx <- res$n == n & res$theta == theta & res$method == m
      res$coverage[idx] <- mean(lo <= theta & theta <= up)
      res$width[idx]    <- mean(up - lo)
    }
  }
}

# ---- tables ----------------------------------------------------------------
wide <- function(var) {
  d <- reshape(res[, c("n", "theta", "method", var)],
               idvar = c("theta", "n"), timevar = "method", direction = "wide")
  names(d) <- sub(paste0(var, "."), "", names(d), fixed = TRUE)
  d <- d[order(d$theta, d$n), ]
  rownames(d) <- NULL
  d
}
tab_cov   <- wide("coverage")
tab_width <- wide("width")
print(tab_cov,   digits = 4)
print(tab_width, digits = 3)

# ---- plot: coverage (solid = no CC, dashed = with CC) ---------------------
cols <- c("red", "red", "blue", "blue", "darkgreen", "darkgreen")
ltys <- c(1, 2, 1, 2, 1, 2)
pchs <- c(16, 1, 17, 2, 15, 0)

plot_res <- function(var, ylim, ylab, hline = NULL) {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for (theta in thetas) {
    plot(NA, xlim = range(ns), ylim = ylim, xaxt = "n",
         xlab = "n", ylab = ylab, main = bquote(theta == .(theta)))
    axis(1, at = ns)
    if (!is.null(hline)) abline(h = hline, lty = 3, col = "gray40")
    for (k in seq_along(methods)) {
      d <- res[res$theta == theta & res$method == names(methods)[k], ]
      lines(d$n, d[[var]], type = "b", col = cols[k], lty = ltys[k],
            pch = pchs[k], lwd = 2)
    }
  }
  plot.new()   # 4th panel: legend
  legend("center", legend = names(methods), col = cols, lty = ltys,
         pch = pchs, lwd = 2, bty = "n", title = "Method")
}

plot_res("coverage", ylim = c(0.6, 1), ylab = "Empirical coverage",
         hline = level)
plot_res("width", ylim = c(0, max(res$width)), ylab = "Mean interval width")
