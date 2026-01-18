if (!interactive()) pdf(NULL)

test_that("balnet works as expected", {
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
    rowMeans(abs(stats.att$control$smd[, -1])),
    pth.att$control$`Mean |SMD|`
  )

  expect_equal(
    apply(abs(stats.att$control$smd[, -1]), 1, max),
    pth.att$control$Lambda,
    tolerance = 1e-4
  )
})
