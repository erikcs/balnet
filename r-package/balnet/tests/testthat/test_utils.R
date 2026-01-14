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

  weights <- runif(n)
  expect_equal(
    col_stats(X, weights, TRUE),
    col_stats(X, weights * 42, TRUE)
  )

  expect_equal(
    col_stats(X, weights)$center[, ],
    apply(X, 2, weighted.mean, weights)
  )

  expect_equal(
    col_stats(X, weights, TRUE)$scale[, ],
    apply(X, 2, function(x) sqrt(sum(weights * (x - weighted.mean(x, weights))^2) / sum(weights)))
  )

  weights <- cbind(1, weights, 1)
  expect_equal(
    col_stats(X, weights, TRUE)$center[1, ],
    col_stats(X, weights, TRUE)$center[3, ]
  )
  expect_equal(
    col_stats(X, weights, TRUE)$scale[1, ],
    col_stats(X, weights, TRUE)$scale[3, ]
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
