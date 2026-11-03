N <- 1e5

par(mfrow = c(1, 1), mar = c(3, 3, 2, 1))

ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)

Y_MAX = c(0.0, 0.0, 0.0, 0.0)

h <- matrix(
  vector("list", length(ns) * length(thetas)),
  nrow=length(ns), ncol=length(thetas)
)

n = ns[4]
p = thetas[2]
S = rbinom(n = 10000, size = n, prob = p)
theta = p * n
theta_hat = mean(S)
p_hat = theta_hat / n
Z = (S/n - p) / sqrt(p_hat * (1 - p_hat) / n)

hist(Z, breaks=20)
