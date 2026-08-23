if (!interactive()) pdf(NULL)

test_that("basic cv.boot.balnet runs", {
  n <- 100
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- cv.boot.balnet(X, W)
  capture.output(print(fit))
  capture.output(summary(fit))
  plot(fit)
  cf <- coef(fit)
  capture.output(print(cf))
  capture.output(summary(cf))
  predict(fit, X)
  wts <- balweights(fit)
  capture.output(print(wts))
  capture.output(summary(wts))

  expect_true(TRUE)
})

test_that("cv.boot.balnet is invariant to sample.weights scale", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  wts <- runif(n)
  seed <- sample.int(1e6, 1)

  set.seed(seed); cv.fit <- cv.boot.balnet(X, W, sample.weights = wts)
  set.seed(seed); cv.fit.scaled <- cv.boot.balnet(X, W, sample.weights = wts * 42)
  expect_equal(
    predict(cv.fit, X),
    predict(cv.fit.scaled, X)
  )
})

test_that("cv.boot.balnet is invariant to W swap", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  seed <- sample.int(1e6, 1)

  set.seed(seed); cv.fit <- cv.boot.balnet(X, W)
  set.seed(seed); cv.fit.swap <- cv.boot.balnet(X, 1 - W)
  expect_equal(
    cv.fit$lambda.min$control,
    cv.fit.swap$lambda.min$treated
  )
  expect_equal(
    cv.fit$lambda.min$treated,
    cv.fit.swap$lambda.min$control
  )
})

test_that("cv.boot.balnet is affine-invariant in X", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  a <- runif(p, 0.5, 10)
  b <- runif(p, -5,  5)
  X.aff <- sweep(sweep(X, 2L, a, `*`), 2L, b, `+`)
  seed <- sample.int(1e6, 1)

  set.seed(seed); cv.fit <- cv.boot.balnet(X, W)
  set.seed(seed); cv.fit.aff <- cv.boot.balnet(X.aff, W)
  expect_equal(
    predict(cv.fit, X),
    predict(cv.fit.aff, X.aff)
  )
})
