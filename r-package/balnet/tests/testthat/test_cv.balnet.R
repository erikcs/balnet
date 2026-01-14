if (!interactive()) pdf(NULL)

test_that("cv.balnet works as expected", {
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

  cv.fit <- cv.balnet(X, W)
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
})
