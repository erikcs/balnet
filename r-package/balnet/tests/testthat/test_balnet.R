if (!interactive()) pdf(NULL)

test_that("basic balnet runs", {
  n <- 100
  p <- 210
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- balnet(X, W)
  capture.output(print(fit))
  plot(fit)
  plot(fit, lambda = 0)
  coef(fit)
  coef(fit, lambda = list(0, 1))
  predict(fit, X)

  fit.gr <- balnet(X, W, groups = list(age = 10:15, 3:7))

  expect_true(TRUE)
})

test_that("balnet is internally consistent (SMD/dev/lmbda)", {
  n <- 111
  p <- 11
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- balnet(X, W)
  capture.output(pth <- print(fit))
  stats <- plot(fit)
  expect_equal(
    rowMeans(abs(stats$control$smd[, -1])),
    pth$control$`Mean |SMD|`
  )
  expect_equal(
    rowMeans(abs(stats$treated$smd[, -1])),
    pth$treated$`Mean |SMD|`
  )
  expect_equal(
    apply(abs(stats$control$smd[, -1]), 1, max),
    pth$control$Lambda,
    tolerance = 1e-4
  )
  expect_equal(
    apply(abs(stats$treated$smd[, -1]), 1, max),
    pth$treated$Lambda,
    tolerance = 1e-4
  )

  fit.att <- balnet(X, W, target = "ATT")
  capture.output(pth.att <- print(fit.att))
  stats.att <- plot(fit.att)
  expect_equal(
    rowMeans(abs(stats.att$smd[, -1])),
    pth.att$`Mean |SMD|`
  )
  expect_equal(
    apply(abs(stats.att$smd[, -1]), 1, max),
    pth.att$Lambda,
    tolerance = 1e-4
  )
})

test_that("balnet is internally consistent (predict/coef)", {
  n <- 111
  p <- 11
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- balnet(X, W)
  expect_equal(
    predict(fit, X),
    predict(fit, X, lambda = fit$lambda)
  )
  expect_equal(
    coef(fit),
    coef(fit, lambda = fit$lambda)
  )
})

test_that("balnet is internally consistent (fits)", {
  n <- 111
  p <- 11
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)

  fit <- balnet(X, W)
  expect_equal(
    predict(fit, X)$control,
    predict(balnet(X, W, target = "control"), X)
  )
  expect_lt(
    mean(abs(predict(fit, X)$control - predict(balnet(X, W, target = "ATT"), X))),
    0.009
  )
  expect_equal(
    predict(fit, X)$treated,
    predict(balnet(X, W, target = "treated"), X)
  )
})

test_that("sample.weighted balnet identical to duplication", {
  n <- 100
  p <- 21
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  to.duplicate <- sample(1:n, 25)
  XX <- rbind(X, X[to.duplicate, ])
  WW <- c(W, W[to.duplicate])
  sample.weights <- rep(1, n)
  sample.weights[to.duplicate] <- 2

  fit.wt <- balnet(X, W, sample.weights = sample.weights)
  fit.dupe <- balnet(XX, WW)

  expect_equal(
    coef(fit.wt),
    coef(fit.dupe)
  )
  expect_equal(
    predict(fit.wt, X),
    predict(fit.dupe, X)
  )
  expect_equal(
    plot(fit.wt),
    plot(fit.dupe)
  )
  expect_equal(
    capture.output(print(fit.wt))[-1],
    capture.output(print(fit.dupe))[-1]
  )
})
