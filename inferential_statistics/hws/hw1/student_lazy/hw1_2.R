# Problem 1.2
# Comparing methods with and without continuity correction (CC)
# Using prop.test and binom.test directly

set.seed(47)

ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)
N <- 1e5
alpha <- 0.10
z <- qnorm(1 - alpha/2)

# Mid-p function for one value of s
get_midp <- function(s, n) {
  lo <- 0; up <- 1
  if (s > 0) {
    f1 <- function(th) 1 - pbinom(s, n, th) + 0.5 * dbinom(s, n, th) - alpha/2
    lo <- uniroot(f1, c(0, 1))$root
  }
  if (s < n) {
    f2 <- function(th) pbinom(s - 1, n, th) + 0.5 * dbinom(s, n, th) - alpha/2
    up <- uniroot(f2, c(0, 1))$root
  }
  c(lo, up)
}

methods <- c("Wald", "Wald+CC", "Score", "Score+CC", "Mid-p", "Exact (C-P)")
cov_res <- data.frame()
w_res   <- data.frame()

for (n in ns) {
  # Compute intervals for all possible s in 0:n
  # 1. Wald
  ci_wald <- sapply(0:n, function(s) {
    p <- s / n
    p + c(-1, 1) * z * sqrt(p * (1 - p) / n)
  })
  
  # 2. Wald + CC
  ci_wald_cc <- sapply(0:n, function(s) {
    p <- s / n
    se <- sqrt(p * (1 - p) / n)
    c(max(0, p - (z * se + 1 / (2 * n))),
      min(1, p + (z * se + 1 / (2 * n))))
  })
  
  # 3. Score (no CC) using prop.test directly
  ci_score <- sapply(0:n, function(s) {
    prop.test(s, n, conf.level = 0.90, correct = FALSE)$conf.int
  })
  
  # 4. Score + CC using prop.test directly
  ci_score_cc <- sapply(0:n, function(s) {
    prop.test(s, n, conf.level = 0.90, correct = TRUE)$conf.int
  })
  
  # 5. Mid-p
  ci_midp <- sapply(0:n, get_midp, n = n)
  
  # 6. Exact Clopper-Pearson using binom.test directly
  ci_exact <- sapply(0:n, function(s) {
    binom.test(s, n, conf.level = 0.90)$conf.int
  })
  
  cis <- list(
    "Wald"        = ci_wald,
    "Wald+CC"     = ci_wald_cc,
    "Score"       = ci_score,
    "Score+CC"    = ci_score_cc,
    "Mid-p"       = ci_midp,
    "Exact (C-P)" = ci_exact
  )
  
  for (theta in thetas) {
    S <- rbinom(N, size = n, prob = theta)
    
    row_c <- data.frame(n = n, theta = theta)
    row_w <- data.frame(n = n, theta = theta)
    
    for (m in methods) {
      lims <- cis[[m]]
      # Column S + 1 corresponds to outcome S
      lo <- lims[1, S + 1]
      up <- lims[2, S + 1]
      row_c[[m]] <- mean(lo <= theta & theta <= up)
      row_w[[m]] <- mean(up - lo)
    }
    cov_res <- rbind(cov_res, row_c)
    w_res   <- rbind(w_res, row_w)
  }
}

cat("Coverage:\n")
print(round(cov_res, 4))
cat("\nMean Width:\n")
print(round(w_res, 3))

# Plots
cols <- c("red", "red", "blue", "blue", "darkgreen", "darkgreen")
ltys <- c(1, 2, 1, 2, 1, 2)
pchs <- c(16, 1, 17, 2, 15, 0)

# Coverage plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (th in thetas) {
  sub_c <- cov_res[cov_res$theta == th, ]
  matplot(sub_c$n, sub_c[, methods], type = "b", pch = pchs, col = cols,
          lty = ltys, lwd = 2, ylim = c(0.6, 1.0), xlab = "n", ylab = "Coverage",
          main = paste("Coverage, theta =", th))
  abline(h = 0.90, lty = 3)
}
plot.new()
legend("center", methods, col = cols, lty = ltys, pch = pchs, lwd = 2, bty = "n")

# Width plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (th in thetas) {
  sub_w <- w_res[w_res$theta == th, ]
  matplot(sub_w$n, sub_w[, methods], type = "b", pch = pchs, col = cols,
          lty = ltys, lwd = 2, ylim = c(0, 0.6), xlab = "n", ylab = "Mean width",
          main = paste("Width, theta =", th))
}
plot.new()
legend("center", methods, col = cols, lty = ltys, pch = pchs, lwd = 2, bty = "n")
