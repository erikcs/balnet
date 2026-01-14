test_that("balnet.fit works", {
  n <- 100
  p <- 210
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)

  fit <- balnet.fit(standardize(X), y)
  predict(fit, X)
  coef(fit)

  expect_true(TRUE)
})

test_that("balnet.fit errs with infeasible lambda", {
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- ifelse(X[, 1] > 0, 1, 0)

  expect_error(
    balnet.fit(standardize(X), y, lambda = 0)
  )
})
