if (!interactive()) pdf(NULL)

test_that("basic cv.balnet runs", {
  n <- 100
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- cv.balnet(X, W)
  capture.output(print(fit))
  plot(fit)
  coef(fit)
  predict(fit, X)

  expect_true(TRUE)
})

test_that("cv.balnet is internally consistent", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  foldid <- sample(rep(seq(10), length.out = nrow(X)))
  cv.fit <- cv.balnet(X, W, foldid = foldid)
  fit <- balnet(X, W)

  expect_equal(
    coef(cv.fit),
    coef(fit, lambda = cv.fit$cv.info$lambda.min)
  )
  expect_equal(
    coef(cv.fit, lambda = list(0, 42)),
    coef(fit, lambda = list(0, 42))
  )
  expect_equal(
    predict(cv.fit, X),
    predict(fit, X, lambda = cv.fit$cv.info$lambda.min)
  )
  expect_equal(
    predict(cv.fit, X, lambda = list(0, 42)),
    predict(fit, X, lambda = list(0, 42))
  )

  expect_equal(
    predict(cv.fit, X)$control,
    predict(cv.balnet(X, W, foldid = foldid, target = "control"), X)
  )
  expect_equal(
    predict(cv.fit, X)$treated,
    predict(cv.balnet(X, W, foldid = foldid, target = "treated"), X)
  )
})

test_that("cv.balnet is invariant to sample.weights scale", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  wts <- runif(n)

  foldid <- sample(rep(seq(3), length.out = nrow(X)))
  cv.fit <- cv.balnet(X, W, foldid = foldid, sample.weights = wts)
  cv.fit.scaled <- cv.balnet(X, W, foldid = foldid, sample.weights = wts * 42)
  expect_equal(
    predict(cv.fit, X),
    predict(cv.fit.scaled, X)
  )
})
