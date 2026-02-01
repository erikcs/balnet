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
    list(control = list(intercepts = 0.680645954386779, betas = new("dgRMatrix",
      p = c(0L, 3L), j = 0:2, Dim = c(1L, 3L), Dimnames = list(
          NULL, NULL), x = c(0.0964925117386255, 0.210324237163201,
      0.145107291953369), factors = list())), treated = list(intercepts = -0.681102012267433,
      betas = new("dgRMatrix", p = c(0L, 3L), j = 0:2, Dim = c(1L,
      3L), Dimnames = list(NULL, NULL), x = c(-0.0496225138747977,
      -0.149449789859248, -0.237897822417655), factors = list())))
  )

  expect_equal(
    coef(fit.att),
    list(intercepts = 0.652873281422005, betas = new("dgRMatrix",
    p = c(0L, 0L), j = integer(0), Dim = c(1L, 3L), Dimnames = list(
        NULL, NULL), x = numeric(0), factors = list()))
  )

})
