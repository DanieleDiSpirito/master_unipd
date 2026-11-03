# HW1 - Problem 1.2
# same as 1.1 but with the continuity correction (CC), see Newcombe (1998)
# - Wald + CC:  p -/+ (z*se + 1/(2n))
# - score + CC: prop.test(correct = TRUE), i.e. Wilson computed at p -/+ 1/(2n)
# - exact: Clopper-Pearson is already the "corrected" one, the uncorrected
#   version is the mid-p interval

set.seed(47)
ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5
alpha <- 0.10
z <- qnorm(1 - alpha/2)

# mid-p interval for a single s (found numerically with uniroot)
midp <- function(s, n) {
  lo <- 0
  up <- 1
  if (s > 0) {
    f <- function(t) 1 - pbinom(s, n, t) + 0.5*dbinom(s, n, t) - alpha/2
    lo <- uniroot(f, c(0, 1))$root
  }
  if (s < n) {
    f <- function(t) pbinom(s - 1, n, t) + 0.5*dbinom(s, n, t) - alpha/2
    up <- uniroot(f, c(0, 1))$root
  }
  c(lo, up)
}

# Wilson interval, used for score and score + CC
wilson <- function(p, n) {
  lo <- (p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
  up <- (p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
  cbind(lo, up)
}

meth <- c("Wald", "Wald+CC", "score", "score+CC", "mid-p", "exact (C-P)")
cover <- array(NA, dim = c(length(ns), length(thetas), length(meth)),
               dimnames = list(ns, thetas, meth))
width <- cover

for (i in 1:length(ns)) {
  n <- ns[i]
  # mid-p for all possible values 0..n (row s+1 = interval for s)
  mp <- t(sapply(0:n, midp, n = n))
  for (j in 1:length(thetas)) {
    theta <- thetas[j]
    S <- rbinom(N, n, theta)
    p <- S/n
    se <- sqrt(p*(1-p)/n)

    lo <- list(); up <- list()
    lo[["Wald"]] <- p - z*se
    up[["Wald"]] <- p + z*se
    lo[["Wald+CC"]] <- pmax(0, p - z*se - 1/(2*n))
    up[["Wald+CC"]] <- pmin(1, p + z*se + 1/(2*n))
    w <- wilson(p, n)
    lo[["score"]] <- w[, 1]
    up[["score"]] <- w[, 2]
    # CC: lower limit from p - 1/(2n), upper limit from p + 1/(2n)
    # (no correction when S = n/2, like prop.test does)
    cc <- ifelse(S == n/2, 0, 0.5/n)
    lo[["score+CC"]] <- ifelse(p - cc <= 0, 0, wilson(p - cc, n)[, 1])
    up[["score+CC"]] <- ifelse(p + cc >= 1, 1, wilson(p + cc, n)[, 2])
    lo[["mid-p"]] <- mp[S + 1, 1]
    up[["mid-p"]] <- mp[S + 1, 2]
    lo[["exact (C-P)"]] <- qbeta(alpha/2, S, n - S + 1)
    up[["exact (C-P)"]] <- qbeta(1 - alpha/2, S + 1, n - S)

    for (m in meth) {
      cover[i, j, m] <- mean(lo[[m]] <= theta & theta <= up[[m]])
      width[i, j, m] <- mean(up[[m]] - lo[[m]])
    }
  }
}

# tables: one for each theta
for (j in 1:length(thetas)) {
  cat("\ntheta =", thetas[j], "\n")
  print(round(cover[, j, ], 4))
}
for (j in 1:length(thetas)) {
  cat("\nmean width, theta =", thetas[j], "\n")
  print(round(width[, j, ], 3))
}

# plots (solid = no CC, dashed = CC)
cols <- c("red", "red", "blue", "blue", "darkgreen", "darkgreen")
ltys <- c(1, 2, 1, 2, 1, 2)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in 1:length(thetas)) {
  matplot(ns, cover[, j, ], type = "b", pch = 1:6, col = cols, lty = ltys,
          lwd = 2, ylim = c(0.6, 1), xlab = "n", ylab = "coverage",
          main = paste("theta =", thetas[j]))
  abline(h = 0.9, lty = 3)
}
plot.new()
legend("center", meth, col = cols, lty = ltys, pch = 1:6, lwd = 2, bty = "n")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (j in 1:length(thetas)) {
  matplot(ns, width[, j, ], type = "b", pch = 1:6, col = cols, lty = ltys,
          lwd = 2, ylim = c(0, 0.6), xlab = "n", ylab = "mean width",
          main = paste("theta =", thetas[j]))
}
plot.new()
legend("center", meth, col = cols, lty = ltys, pch = 1:6, lwd = 2, bty = "n")
