EXPECTED_AD_ERROR <- "adelie_core solver: max coordinate descents reached at lambda index: 0."

#' Low-level fit function for adelie cbps solver.
#'
#' @param stan List containing the standardized feature matrix along with mean and scales.
#' @param y The 0/1 outcome.
#' @param weights Sample weights.
#' @param target_scale Gradient scaling for glm.
#' @param lambda Optional `lambda` sequence. By default, the program computes `lambda` sequence based on `lmda_path_size` and `min_ratio`.
#' @param lmda_path_size Number of values for `lambda`, if generated automatically. Default is 100.
#' @param min_ratio Ratio between smallest and largest value of lambda. Default is 1e-2.
#' @param penalty Penalty factor per feature. Default is 1.
#' @param groups List of group indices. Default is each variable is a group.
#' @param alpha Elastic net mixing parameter. Default is 1 (lasso). 0 is ridge.
#' @param irls_max_iters Maximum number of IRLS iterations, default is 1e4.
#' @param irls_tol IRLS convergence tolerance, default is 1e-7.
#' @param max_iters Maximum total number of coordinate descent iterations at each BASIL step, default is 1e5.
#' @param tol Coordinate descent convergence tolerance, default 1e-7.
#' @param newton_max_iters Maximum number of iterations for the BCD update, default 1000.
#' @param newton_tol Convergence tolerance for the BCD update, default 1e-12.
#' @param screen_rule Screen rule, with default `"pivot"`.
#' @param max_screen_size Maximum number of screen groups.
#' @param max_active_size Maximum number of active groups.
#' @param pivot_subset_ratio Subset ratio of pivot rule.
#' @param pivot_subset_min Minimum subset of pivot rule.
#' @param pivot_slack_ratio Slack ratio of pivot rule.
#' @param progress_bar Progress bar. Default is `FALSE`.
#' @param n_threads Number of threads, default 1.
#'
#' @return A balnet.fit object.
#'
#' @keywords internal
#' @export
balnet.fit <- function(
  stan,
  y,
  weights = NULL,
  target_scale = 1,
  lambda = NULL,
  lmda_path_size = 100L,
  min_ratio = 1e-2,
  penalty = NULL,
  groups = NULL,
  alpha = 1.0,
  irls_max_iters = as.integer(1e4),
  irls_tol = 1e-7,
  max_iters = as.integer(1e5),
  tol = 1e-7,
  newton_max_iters = 1000L,
  newton_tol = 1e-12,
  screen_rule = c("pivot", "strong"),
  max_screen_size = NULL,
  max_active_size = NULL,
  pivot_subset_ratio = 0.1,
  pivot_subset_min = 1L,
  pivot_slack_ratio = 1.25,
  progress_bar = FALSE,
  n_threads = 1L
)
{
  screen_rule <- match.arg(screen_rule)
  X <- stan[["X"]]
  groups <- process_groups(groups, ncol(X))
  if (!is.numeric(lmda_path_size) || lmda_path_size <= 0) {
    stop("Invalid lambda path size.")
  }
  if (!is.numeric(min_ratio) || min_ratio <= 0 || min_ratio >= 1) {
    stop("Invalid lambda ratio.")
  }
  if (!is.null(lambda) && (!is.numeric(lambda) || length(lambda) == 0)) {
    stop("Invalid 'lambda'.")
  }
  if (!is.null(penalty) && (!is.numeric(penalty) || length(penalty) != groups$G) || any(penalty < 0)) {
    stop("Invalid penalty factor.")
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("Invalid 'alpha'.")
  }
  if (!is.numeric(n_threads) || n_threads <= 0) {
    stop("Invalid number of threads.")
  }
  if (
    any(!is.finite(c(irls_max_iters, irls_tol, max_iters, tol, newton_max_iters, newton_tol))) ||
    any(c(irls_max_iters, irls_tol, max_iters, tol, newton_max_iters, newton_tol) < 0)
  ) {
    stop("Invalid tolerance/iteration input.")
  }

  if (is.null(lambda)) {
    lmda_path <- double(0)
    setup_lmda_path <- TRUE
  } else {
    lmda_path <- sort(lambda, decreasing = TRUE)
    setup_lmda_path <- FALSE
    lmda_path_size <- length(lmda_path)
  }
  lmda_max <- -1.0 # Let the solver compute this, is simply max(abs(crossprod(X, weights * (y / mean(y) - 1)))) (or of colMeans(X * y)), but does it faster
  setup_lmda_max <- TRUE

  if (is.null(penalty)) {
    penalty <- sqrt(groups$group_sizes)
  }
  if (is.null(max_screen_size)) {
    max_screen_size <- groups$G
  }
  if (is.null(max_active_size)) {
    max_active_size <- groups$G
  }
  max_screen_size <- min(max_screen_size, groups$G)
  max_active_size <- min(max_active_size, groups$G)

  # GLM args
  if (is.null(weights)) {
    weights <- rep_len(1 / length(y), length(y))
  } else {
    weights <- weights / sum(weights)
  }
  y <- as.double(y)

  beta0 <- 0.0
  offsets <- double(length(y)) # unused
  eta <- as.double(offsets)
  resid <- weights * (y * exp(-eta) - (1 - y)) # negative gradient of calib. loss
  grad <- double(ncol(X)) # This is crossprod(X, resid), but initialized in make_state since that's faster

  # *Hardcoded solver defaults*
  intercept <- TRUE

  # Disable deviance based CD tolerance by defining loss_null = 1 and loss_full = 0
  # (loss_null/loss_full is not used anywhere else in the solver, for path progress printout, deviance is re-defined in C++)
  loss_null <- 1.0
  loss_full <- 0.0
  setup_loss_null <- FALSE
  # Disable deviance-based early exit (not relevant for cbps)
  adev_tol <- 0.9
  ddev_tol <- 0
  early_exit <- FALSE

  # *Unused solver options*
  constraints <- replicate(groups$G, NULL, FALSE)
  dual_groups <- rep(0L, groups$G)

  # Unused warm start options
  lmda <- Inf
  screen_set <- (0:(groups$G-1))[(penalty <= 0) | (alpha <= 0)]
  screen_beta <- double(sum(groups$group_sizes[screen_set + 1]))
  screen_is_active <- as.integer(rep_len(1, length(screen_set)))
  active_set_size <- length(screen_set)
  active_set <- integer(groups$G)
  if (active_set_size > 0) {
    active_set[1:active_set_size] <- 0:(active_set_size-1)
  }

  args <- list(
    # State args
    "X" = X,
    "eta" = eta,
    "resid" = resid,
    "constraints" = constraints,
    "groups" = groups$groups,
    "group_sizes" = groups$group_sizes,
    "dual_groups" = dual_groups,
    "alpha" = alpha,
    "penalty" = penalty,
    "offsets" = offsets,
    "lmda_path" = lmda_path,
    "loss_null" = loss_null,
    "loss_full" = loss_full,
    "lmda_max" = lmda_max,
    "min_ratio" = min_ratio,
    "lmda_path_size" = lmda_path_size,
    "max_screen_size" = max_screen_size,
    "max_active_size" = max_active_size,
    "pivot_subset_ratio" = pivot_subset_ratio,
    "pivot_subset_min" = pivot_subset_min,
    "pivot_slack_ratio" = pivot_slack_ratio,
    "screen_rule" = screen_rule,
    "irls_max_iters" = irls_max_iters,
    "irls_tol" = irls_tol,
    "max_iters" = max_iters,
    "tol" = tol,
    "adev_tol" = adev_tol,
    "ddev_tol" = ddev_tol,
    "newton_tol" = newton_tol,
    "newton_max_iters" = newton_max_iters,
    "early_exit" = early_exit,
    "setup_loss_null" = setup_loss_null,
    "setup_lmda_max" = setup_lmda_max,
    "setup_lmda_path" = setup_lmda_path,
    "intercept" = intercept,
    "n_threads" = n_threads,
    "screen_set" = screen_set,
    "screen_beta" = screen_beta,
    "screen_is_active" = screen_is_active,
    "active_set_size" = active_set_size,
    "active_set" = active_set,
    "beta0" = beta0,
    "lmda" = lmda,
    "grad" = grad,
    # GLM args
    "y" = y,
    "weights" = weights,
    "target_scale" = target_scale,
    # Solver args
    "progress_bar" = progress_bar
  )
  fit <- rcpp_solver(args)
  class(fit) <- "balnet.fit"

  fit[["lmda_path"]] <- drop(fit[["lmda_path"]])
  fit[["stan"]] <- stan[-1]

  error <- fit[["error"]]
  if (nzchar(error)) {
    if (!is.null(lambda)) {
      if (length(fit[["lmdas"]]) == 0) {
        stop(paste("None of the supplied lambdas converged.",
          "Try increasing the magnitude or let the program calculate the lambda path automatically."))
      }
    } else {
      # Expected exception when calibration loss starts to diverge at end of lambda path.
      # Acts as a natural stopping criterion, returning the path up to but not including this lambda.
      if (!identical(error, EXPECTED_AD_ERROR)) {
        warning(error)
      }
    }
  }

  fit
}

#' Extract coefficients from a balnet.fit object.
#'
#' @param object A `balnet.fit` object.
#' @param lambda Value(s) for the penalty parameter. If NULL (default), the
#'   lambda path on which `object` was fit is used. If different lambda values
#'   are supplied, linear interpolation is used, as in `glmnet`. Note: no
#'   extrapolation is performed.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Coefficients.
#'
#' @keywords internal
#' @export
coef.balnet.fit <- function(
  object,
  lambda = NULL,
  ...
)
{
  intercepts <- object[["intercepts"]]
  betas <- object[["betas"]]
  if (is.null(lambda)) {
    lambda <- object[["lmdas"]]
  } else {
    lamlist <- lambda.interp(object[["lmdas"]], lambda)
    betas <- Matrix::Diagonal(x = lamlist$frac) %*% betas[lamlist$left, , drop = FALSE] +
      Matrix::Diagonal(x = 1 - lamlist$frac) %*% betas[lamlist$right, , drop = FALSE]
    intercepts <- diag(x = lamlist$frac, nrow = length(lambda)) %*% intercepts[lamlist$left] +
      diag(x = 1 - lamlist$frac, nrow = length(lambda)) %*% intercepts[lamlist$right]
  }

  # Coefficients kept on unit scale and re-scaled when needed
  stan <- object[["stan"]]
  intercepts <- intercepts - betas %*% (stan[["center"]] / stan[["scale"]])
  betas <- betas %*% Matrix::Diagonal(x = 1 / stan[["scale"]])

  list(intercepts = drop(intercepts), betas = betas)
}

#' Predict using a balnet.fit object.
#'
#' @param object A balnet.fit object.
#' @param newx A numeric matrix.
#' @param lambda Value(s) for the penalty parameter. If NULL (default), the
#'   lambda path on which `object` was fit is used. If different lambda values
#'   are supplied, linear interpolation is used, as in `glmnet`. Note: no
#'   extrapolation is performed.
#' @param type The type of predictions.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Predictions.
#'
#' @keywords internal
#' @export
predict.balnet.fit <- function(
  object,
  newx,
  lambda = NULL,
  type = c("response", "link"),
  ...
)
{
  type <- match.arg(type)

  coefs <- coef(object, lambda = lambda)
  intercepts <- coefs[["intercepts"]]
  betas <- coefs[["betas"]]
  eta <- tcrossprod(newx, betas) + matrix(intercepts, nrow(newx), length(intercepts), byrow = TRUE)

  if (type == "response") {
    out <- 1 / (1 + exp(-eta))
  } else if (type == "link") {
    out <- eta
  }

  out
}
