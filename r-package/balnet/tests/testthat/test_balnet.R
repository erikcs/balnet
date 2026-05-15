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
  W <- rbinom(n, 1, 0.3)

  fit <- balnet(X, W)
  capture.output(pth <- print(fit))
  stats.pth <- plot(fit)
  stats.smd <- plot(fit, lambda = 0)

  # control
  start.end0 <- c(1, length(pth$control$`Avg|SMD|`))
  expect_equal(
    unname(colMeans(abs(stats.smd$control))),
    pth$control$`Avg|SMD|`[start.end0],
  )
  expect_equal(
    unname(apply(abs(stats.smd$control), 2, max)),
    pth$control$Lambda[start.end0],
    tolerance = 2e-4
  )
  expect_equal(
    unname((1 - colSums(abs(stats.smd$control)) / sum(abs(stats.smd$control$lambda.max))) * 100),
    stats.pth$control$pbr[start.end0]
  )

  # treated
  start.end1 <- c(1, length(pth$treated$`Avg|SMD|`))
  expect_equal(
    unname(colMeans(abs(stats.smd$treated))),
    pth$treated$`Avg|SMD|`[start.end1],
  )
  expect_equal(
    unname(apply(abs(stats.smd$treated), 2, max)),
    pth$treated$Lambda[start.end1],
    tolerance = 2e-4
  )
  expect_equal(
    unname((1 - colSums(abs(stats.smd$treated)) / sum(abs(stats.smd$treated$lambda.max))) * 100),
    stats.pth$treated$pbr[start.end1]
  )

  fit.att <- balnet(X, W, target = "ATT")
  capture.output(pth.att <- print(fit.att))
  stats.pth.att <- plot(fit.att)
  stats.smd.att <- plot(fit.att, lambda = 0)

  start.end <- c(1, length(pth.att$`Avg|SMD|`))
  expect_equal(
    unname(colMeans(abs(stats.smd.att))),
    pth.att$`Avg|SMD|`[start.end],
  )
  expect_equal(
    unname(apply(abs(stats.smd.att), 2, max)),
    pth.att$Lambda[start.end0],
    tolerance = 2e-4
  )
  expect_equal(
    unname((1 - colSums(abs(stats.smd.att)) / sum(abs(stats.smd.att$lambda.max))) * 100),
    stats.pth.att$pbr[start.end]
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

  Y <- runif(n)
  YY <- c(Y, Y[to.duplicate])
  expect_equal(
    apply(balweights(fit.wt)$control, 2, function(x) weighted.mean(Y, x)),
    apply(balweights(fit.dupe)$control, 2, function(x) weighted.mean(YY, x))
  )
  expect_equal(
    apply(balweights(fit.wt)$treated, 2, function(x) weighted.mean(Y, x)),
    apply(balweights(fit.dupe)$treated, 2, function(x) weighted.mean(YY, x))
  )
})

test_that("balnet has not changed", {
  set.seed(42)
  n <- 111
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.3)

  fit <- balnet(X, W, nlambda = 5)
  fit.att <- balnet(X, W, target = "ATT", nlambda = 5)

  expect_equal(
    coef(fit, lambda = fit[["_lambda"]]),
    list(control = new("dgCMatrix", i = c(0L, 0L, 1L, 2L, 3L, 0L,
    1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L), p = c(0L, 1L, 5L,
    9L, 13L, 17L), Dim = 4:5, Dimnames = list(c("(Intercept)", "X1",
    "X2", "X3"), NULL), x = c(0.652873281422005, 0.66802481765921,
    0.0507573767431381, 0.158769917286064, 0.103686341052803, 0.680645954386779,
    0.0964925117386255, 0.210324237163201, 0.145107291953369, 0.685941154728131,
    0.111004789201444, 0.227049604256553, 0.157882123136175, 0.687747421060082,
    0.115601582272039, 0.232383701326134, 0.161889105216467), factors = list()),
        treated = new("dgCMatrix", i = c(0L, 0L, 1L, 2L, 3L, 0L,
        1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L), p = c(0L, 1L,
        5L, 9L, 13L, 17L), Dim = 4:5, Dimnames = list(c("(Intercept)",
        "X1", "X2", "X3"), NULL), x = c(-0.652873281422005, -0.682941839800975,
        -0.00207812408065897, -0.118494422306579, -0.156693688472859,
        -0.682398475582524, -0.0379564149766521, -0.142215753917237,
        -0.218193176446944, -0.681102012267433, -0.0496225138747977,
        -0.149449789859248, -0.237897822417655, -0.680577640484851,
        -0.0533492030403526, -0.151711346250273, -0.244157630108172
        ), factors = list()))
  )

  expect_equal(
    coef(fit.att, lambda = fit.att[["_lambda"]]),
    new("dgCMatrix", i = c(0L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L,
    1L, 2L, 3L, 0L, 1L, 2L, 3L), p = c(0L, 1L, 5L, 9L, 13L, 17L),
        Dim = 4:5, Dimnames = list(c("(Intercept)", "X1", "X2", "X3"
        ), NULL), x = c(0.652873281422005, 0.667336389848779, 0.0527093774420999,
        0.152995377993314, 0.112631732214514, 0.680472124137364,
        0.0970443032964536, 0.208453179362447, 0.147901280091804,
        0.68588965859546, 0.111172722042331, 0.226453080921526, 0.158762034363933,
        0.687731452647602, 0.115654028489226, 0.232194567332454,
        0.162166992288202), factors = list())
  )

  expect_equal(
    coef(fit.att, lambda = c(-100, fit.att$lambda * 0.942, 0, 42)),
    new("dgCMatrix", i = c(0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L,
    2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L, 2L, 3L, 0L, 1L,
    2L, 3L, 0L), p = c(0L, 4L, 8L, 12L, 16L, 20L, 24L, 28L, 29L),
        Dim = c(4L, 8L), Dimnames = list(c("(Intercept)", "X1", "X2",
        "X3"), NULL), x = c(0.687731452647602, 0.115654028489226,
        0.232194567332454, 0.162166992288202, 0.654100093870745,
        0.00447099741654189, 0.0129776137178316, 0.00955382529996484,
        0.668450609688184, 0.0564700240021002, 0.157699506531369,
        0.115623421076997, 0.680931658735221, 0.0982427261725733,
        0.209979995220917, 0.148822528011033, 0.686045886136906,
        0.111552842440619, 0.226940094279606, 0.159050855041043,
        0.687731452647602, 0.115654028489226, 0.232194567332454,
        0.162166992288202, 0.687731452647602, 0.115654028489226,
        0.232194567332454, 0.162166992288202, 0.652873281422005),
        factors = list())
  )
})
