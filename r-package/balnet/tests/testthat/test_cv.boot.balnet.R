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

  fit <- cv.boot.balnet(X, W, type.measure = "imbalance.inf")

  expect_true(TRUE)
})

test_that("cv.boot.balnet is invariant to sample.weights scale", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  wts <- runif(n)
  seed <- sample.int(1e6, 1)

  for (tm in c("imbalance.mean", "imbalance.inf")) {
    set.seed(seed)
    cv.fit <- cv.boot.balnet(X, W, sample.weights = wts,
                             type.measure = tm, B = 20)

    set.seed(seed)
    cv.fit.scaled <- cv.boot.balnet(X, W, sample.weights = wts * 42,
                                    type.measure = tm, B = 20)

    expect_equal(
      predict(cv.fit, X),
      predict(cv.fit.scaled, X),
      label = paste0("predict(), type.measure = ", tm)
    )
  }
})

test_that("cv.boot.balnet is invariant to W swap", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  seed <- sample.int(1e6, 1)

  for (tm in c("imbalance.mean", "imbalance.inf")) {
    set.seed(seed)
    cv.fit <- cv.boot.balnet(X, W, type.measure = tm, B = 20)

    set.seed(seed)
    cv.fit.swap <- cv.boot.balnet(X, 1 - W, type.measure = tm, B = 20)

    expect_equal(
      cv.fit$lambda.min$control,
      cv.fit.swap$lambda.min$treated,
      label = paste0("lambda.min control->treated, type.measure = ", tm)
    )
    expect_equal(
      cv.fit$lambda.min$treated,
      cv.fit.swap$lambda.min$control,
      label = paste0("lambda.min treated->control, type.measure = ", tm)
    )
  }
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

  for (tm in c("imbalance.mean", "imbalance.inf")) {
    set.seed(seed)
    cv.fit <- cv.boot.balnet(X, W, type.measure = tm, B = 20)

    set.seed(seed)
    cv.fit.aff <- cv.boot.balnet(X.aff, W, type.measure = tm, B = 20)

    expect_equal(
      predict(cv.fit, X),
      predict(cv.fit.aff, X.aff),
      label = paste0("predict(), type.measure = ", tm)
    )
  }
})
