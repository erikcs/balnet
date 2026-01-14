library(balnet)
n <- 100
p <- 30
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
gr <- list(5:15, 17:19)

fit1 <- balnet(X, W, penalty.factor = runif(ncol(X)))

fit2 <- balnet(X, W, groups = gr, alpha = 0.5)
