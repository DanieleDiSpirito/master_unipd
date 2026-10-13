hdd <- read.table('data/hdds.tab', sep='\t', header=TRUE)
colSums(hdd)

plot(
  function(x) dbinom(as.integer(x), 400, prob=0.2),
  xlim=c(0,400*0.2*2),
  type="h",
  ylab="Probability",
  xlab="# used HDDs",
)
plot(
  function(x) dnorm(x, mean=0.2 * 400, sd = sqrt(0.2*0.8*400)),
  xlim=c(0, 400*0.2*2),
  add=TRUE,
  n = 100,
)

S = rbinom(n = 100e5, size=400, prob=0.2)
W = (S - 400*0.2)/sqrt(0.2*0.8*400)
hist(
  W,
  freq=FALSE,
  xlab="Std # used HDDs",
  ylab="Density",
  main="W vs N",
  xlim=c(-5,5),
  breaks=64,
)
plot(
  dnorm,
  xlim=c(-5, 5),
  add=TRUE,
  n=100,
)

