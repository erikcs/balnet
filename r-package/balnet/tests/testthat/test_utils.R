test_that("col stats works as expected", {
  n <- 101
  p <- 13
  X <- matrix(rnorm(n * p), n, p)

  weights <- rep(1, n)
  expect_equal(
    col_stats(X, weights, TRUE),
    col_stats(X, weights / n, TRUE)
  )
  expect_equal(
    col_stats(X, weights, TRUE),
    col_stats(X, weights * 42, TRUE)
  )

  weights.rn <- runif(n)
  expect_equal(
    col_stats(X, weights.rn, TRUE),
    col_stats(X, weights.rn * 42, TRUE)
  )

  expect_equal(
    col_stats(X, weights.rn)$center[, ],
    apply(X, 2, weighted.mean, weights.rn)
  )

  expect_equal(
    col_stats(X, weights.rn, TRUE)$scale[, ],
    apply(X, 2, function(x) sqrt(sum(weights.rn * (x - weighted.mean(x, weights.rn))^2) / sum(weights.rn)))
  )

  weights.mat <- cbind(1, weights.rn, 1, runif(n))
  expect_equal(
    col_stats(X, weights.mat, TRUE)$center[1, ],
    col_stats(X, weights.mat, TRUE)$center[3, ]
  )
  expect_equal(
    col_stats(X, weights.mat, TRUE)$scale[1, ],
    col_stats(X, weights.mat, TRUE)$scale[3, ]
  )

  expect_equal(
    col_stats(X, weights.mat, TRUE)$center[4, ],
    col_stats(X, weights.mat[, 4], TRUE)$center[, ]
  )
  expect_equal(
    col_stats(X, weights.mat, TRUE)$scale[4, ],
    col_stats(X, weights.mat[, 4], TRUE)$scale[1, ]
  )

  n <- 1000
  p <- 10
  X <- matrix(rnorm(n * p), n, p) * 4.42
  W <- rbinom(n, 1, 1 / (1 + exp(2.5 - X[, 1])))
  X[W == 1, 8] <- 1.42

  expect_equal(
    col_stats(X, W, compute_sd = TRUE)$scale[, 8],
    0
  )
})

test_that("standardize works as expected", {
  R_scale <- function(X, stan) {
    sc.x <- scale(X, stan$center, stan$scale)
    attributes(sc.x)$`scaled:center` <- NULL
    attributes(sc.x)$`scaled:scale` <- NULL

    sc.x
  }
  n <- 101
  p <- 13
  X <- matrix(rnorm(n * p), n, p) * 42

  stan <- standardize(X)
  expect_equal(stan$X, R_scale(X, stan))

  X.cpy <- X * 1.0
  stan <- standardize(X.cpy, inplace = TRUE)
  expect_equal(X.cpy, stan$X)
})

test_that("collapse mat works as expected", {
  n <- 101
  p <- 133
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- make.names(1:p)

  expect_equal(
    unname(collapse_X(X, as.list(1:p), colnames(X))),
    unname(X)
  )

  x1 <- collapse_X(X, list(c(100:133), 5:10), colnames(X))
  expect_equal(
    x1[, "Grp_100-133"],
    rowMeans(X[, 100:133])
  )
  expect_equal(
    x1[, "Grp_5-10"],
    rowMeans(X[, 5:10])
  )

  expect_equal(
    collapse_X(X, list(c(100:133), age=5:10), colnames(X))[, "age"],
    rowMeans(X[, 5:10])
  )
})

test_that("sp_tcrossprod_plus works as expected", {
  n <- 50
  p <- 542
  L <- 100
  X <- matrix(rnorm(n * p), n, p)
  beta <- Matrix::rsparsematrix(nrow = L, ncol = p, density = 0.05)
  beta <- as(beta, "RsparseMatrix")
  intercepts <- rnorm(L)

  expect_equal(
    sp_tcrossprod_plus(X, beta, intercepts),
    as.matrix(tcrossprod(X, beta) + matrix(intercepts, nrow(X), length(intercepts), byrow = TRUE))
  )
})
