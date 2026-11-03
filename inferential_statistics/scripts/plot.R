N <- 1e5

ns <- c(10, 30, 50, 100)
thetas <- c(0.1, 0.5, 0.9)

par(mfrow = c(4, 3), mar = c(3, 3, 2, 1))

Y_MAX = c(0.0, 0.0, 0.0, 0.0)

h <- matrix(
  vector("list", length(ns) * length(thetas)),
  nrow=length(ns), ncol=length(thetas)
)

set.seed(47)
for(i in seq_along(ns)) {
  for(j in seq_along(thetas)) {
    S <- rbinom(n = N, size = ns[i], prob = thetas[j])
    h[[i, j]] <- hist(
      S, breaks = seq(min(S), max(S), by = 1), 
      right=TRUE, plot=FALSE
    )
    Y_MAX[i] <- max(max(h[[i, j]]$density), Y_MAX[i]);
  }
}

for(i in seq_along(ns)) {
  for(j in seq_along(thetas)) {
    plot(h[[i,j]],
         xlim = c(0, ns[i]),
         ylim = c(0, Y_MAX[i]),
         freq=FALSE,
         main=bquote(~ n == .(ns[i]) ~ "," ~ theta == .(thetas[j]))
    );
  }
}

