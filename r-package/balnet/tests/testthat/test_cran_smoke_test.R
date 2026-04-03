test_that("basic balnet runs", {
  n <- 100
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- balnet(X, W)

  expect_true(TRUE)
})
