tryCatch(
  {
    attachNamespace("CVXR")
  },
  error = function(e) {
    install.packages("CVXR", repos = "http://cran.us.r-project.org")
    attachNamespace("CVXR")
  }
)

cvx_mu1_lasso = function(y, X, lambda, standardize = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  if (standardize) {
    xm <- colMeans(X)
    xs <- apply(X, 2, sd) * sqrt((n - 1) / n)
    xs[xs == 0] <- 1

    Xs <- scale(X, center = xm, scale = xs)
  } else {
    xm <- rep(0, p)
    xs <- rep(1, p)
    Xs <- X
  }

  beta0 <- Variable(1)
  beta <- Variable(p)
  eta <- beta0 + Xs %*% beta

  loss_terms <- y * exp(-eta) + (1 - y) * eta
  loss <- mean(loss_terms)
  objective <- loss + lambda * p_norm(beta, 1)

  problem <- Problem(Minimize(objective))
  result <- solve(problem)

  beta_std  <- as.numeric(result$getValue(beta))
  beta0_std <- as.numeric(result$getValue(beta0))

  # Convert back to original scale
  beta_orig  <- beta_std / xs
  beta0_orig <- beta0_std - sum(beta_std * xm / xs)

  eta_hat <- as.numeric(beta0_orig + X %*% beta_orig)
  ps_hat  <- 1 / (1 + exp(-eta_hat))

  list(
    beta0 = beta0_orig,
    beta = beta_orig,
    ps = ps_hat,
    eta = eta_hat,
    beta0_std = beta0_std,
    beta_std = beta_std,
    status = result$status,
    value = result$value,
    loss = mean(y * exp(-eta_hat) + (1 - y) * eta_hat)
  )
}

# Pick first, mid, and 4-th to last lmbda in path
# (in some settings, the final one is too low for CVX to converge)
get_example_lambdas = function(fit) {
  v = fit$lmdas
  v[c(1, ceiling(length(v) / 2), max(length(v) - 4, 1))]
}

test_that("works as expected vs CVXR", {
  # data 1: n > p and overlap
  setup = "data1"
  n = 100
  p = 10
  X = 15 * matrix(rnorm(n * p), n, p); X[, 2] = X[, 2] * 42; X[, 3] = rbinom(n, 1, 0.5)
  y = rbinom(n, 1, 0.5)
  fit = balnet.fit(standardize(X), y)
  for (lambda in get_example_lambdas(fit)) {
    test_that(paste("setup =", setup, "lambda =", lambda), {
      pp = predict(fit, newx = X, type = "response", lambda = lambda)[, ]
      fit.cvx = cvx_mu1_lasso(y, X, lambda)

      expect_lt(max(abs(pp - fit.cvx$ps)), 1e-3)
      expect_lt(sqrt(mean((pp - fit.cvx$ps)^2)), 1e-3)
    })
  }

  # data 2: n < p and overlap
  setup = "data2"
  n = 100
  p = 200
  X = 15 * matrix(rnorm(n * p), n, p); X[, 2] = X[, 2] * 42; X[, 3] = rbinom(n, 1, 0.5)
  y = rbinom(n, 1, 0.5)
  fit = balnet.fit(standardize(X), y)
  for (lambda in get_example_lambdas(fit)) {
    test_that(paste("setup =", setup, "lambda =", lambda), {
      pp = predict(fit, newx = X, type = "response", lambda = lambda)[, ]
      fit.cvx = cvx_mu1_lasso(y, X, lambda)

      expect_lt(max(abs(pp - fit.cvx$ps)), 1e-4)
      expect_lt(sqrt(mean((pp - fit.cvx$ps)^2)), 1e-4)
    })
  }

  # data 3: n > p and low overlap
  setup = "data3"
  n = 100
  p = 10
  X = matrix(rnorm(n * p), n, p); X[, 2] = X[, 2] * 42; X[, 3] = rbinom(n, 1, 0.5)
  y = rbinom(n, 1, 1 / (1 + exp(2.5 - X[, 1])))

  fit = balnet.fit(standardize(X), y)
  for (lambda in get_example_lambdas(fit)[1:2]) {
    test_that(paste("setup =", setup, "lambda =", lambda), {
      pp = predict(fit, newx = X, type = "response", lambda = lambda)[, ]
      fit.cvx = cvx_mu1_lasso(y, X, lambda)

      expect_lt(max(abs(pp - fit.cvx$ps)), 1e-4)
      expect_lt(sqrt(mean((pp - fit.cvx$ps)^2)), 1e-4)
    })
  }

  # data 4: n < p and low overlap
  setup = "data4"
  n = 100
  p = 200
  X = matrix(rnorm(n * p), n, p); X[, 2] = X[, 2] * 42; X[, 3] = rbinom(n, 1, 0.5)
  y = rbinom(n, 1, 1 / (1 + exp(2.5 - X[, 1])))

  fit = balnet.fit(standardize(X), y)
  for (lambda in get_example_lambdas(fit)[1:2]) {
    test_that(paste("setup =", setup, "lambda =", lambda), {
      pp = predict(fit, newx = X, type = "response", lambda = lambda)[, ]
      fit.cvx = cvx_mu1_lasso(y, X, lambda)

      expect_lt(max(abs(pp - fit.cvx$ps)), 1e-4)
      expect_lt(sqrt(mean((pp - fit.cvx$ps)^2)), 1e-4)
    })
  }

  # data 5: n < p and complete separation
  setup = "data5"
  n = 100
  p = 200
  X = 15 * matrix(rnorm(n * p), n, p); X[, 2] = X[, 2] * 42; X[, 3] = rbinom(n, 1, 0.5)
  y = ifelse(X[, 1] >= 0, 1, 0)

  fit = balnet.fit(standardize(X), y)
  for (lambda in get_example_lambdas(fit)[1:2]) {
    test_that(paste("setup =", setup, "lambda =", lambda), {
      pp = predict(fit, newx = X, type = "response", lambda = lambda)[, ]
      fit.cvx = cvx_mu1_lasso(y, X, lambda)

      expect_lt(max(abs(pp - fit.cvx$ps)), 1e-4)
      expect_lt(sqrt(mean((pp - fit.cvx$ps)^2)), 1e-4)
    })
  }

  expect_true(TRUE)
})
