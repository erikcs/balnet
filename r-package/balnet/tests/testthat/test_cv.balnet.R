if (!interactive()) pdf(NULL)

test_that("basic cv.balnet runs", {
  n <- 100
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- cv.balnet(X, W)
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

  fit <- cv.balnet(X, W, type.measure = "imbalance.mean")
  fit <- cv.balnet(X, W, type.measure = "imbalance.inf")

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
    coef(fit, lambda = cv.fit$`_cv.info`$lambda.min)
  )
  expect_equal(
    coef(cv.fit, lambda = list(0, 42)),
    coef(fit, lambda = list(0, 42))
  )
  expect_equal(
    predict(cv.fit, X),
    predict(fit, X, lambda = cv.fit$`_cv.info`$lambda.min)
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

test_that("cv.balnet has not changed", {
  set.seed(42)
  n <- 111
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.3)

  fit <- cv.balnet(X, W, nlambda = 5)
  fit.att <- cv.balnet(X, W, target = "ATT", nlambda = 5)

  expect_equal(
    coef(fit),
    structure(list(control = new("dgCMatrix", i = 0L, p = 0:1, Dim = c(4L,
1L), Dimnames = list(c("(Intercept)", "X1", "X2", "X3"), NULL),
    x = 0.652873281422005, factors = list()), treated = new("dgCMatrix",
    i = 0L, p = 0:1, Dim = c(4L, 1L), Dimnames = list(c("(Intercept)",
    "X1", "X2", "X3"), NULL), x = -0.652873281422005, factors = list())), class = "coef.balnet.contrast")
  )

  expect_equal(
    coef(fit.att),
    new("dgCMatrix", i = 0L, p = 0:1, Dim = c(4L, 1L), Dimnames = list(
      c("(Intercept)", "X1", "X2", "X3"), NULL), x = 0.652873281422005,
      factors = list())
  )

})

test_that("cv.balnet is invariant to W swap", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  foldid <- sample(rep(seq(3), length.out = nrow(X)))
  for (crit in c("balance.loss", "imbalance.mean", "imbalance.inf")) {
    cv.fit <- cv.balnet(X, W, foldid = foldid, type.measure = crit)
    cv.fit.swap <- cv.balnet(X, 1 - W, foldid = foldid, type.measure = crit)
    expect_equal(
      cv.fit$lambda.min$control,
      cv.fit.swap$lambda.min$treated
    )
    expect_equal(
      cv.fit$lambda.min$treated,
      cv.fit.swap$lambda.min$control
    )
  }
})

test_that("cv.balnet is affine-invariant in X", {
  n <- 111
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  a <- runif(p, 0.5, 10)
  b <- runif(p, -5,  5)
  X.aff <- sweep(sweep(X, 2L, a, `*`), 2L, b, `+`)
  foldid <- sample(rep(seq(3), length.out = nrow(X)))
  for (crit in c("balance.loss", "imbalance.mean", "imbalance.inf")) {
    cv.fit <- cv.balnet(X, W, foldid = foldid, type.measure = crit)
    cv.fit.aff <- cv.balnet(X.aff, W, foldid = foldid, type.measure = crit)
    expect_equal(
      predict(cv.fit, X),
      predict(cv.fit.aff, X.aff)
    )
  }
})
