#' Pathwise estimation of covariate balancing propensity scores.
#'
#' Fits regularized logistic regression models using covariate balancing loss functions, targeting
#' the ATE, ATT, or treated/control means.
#'
#' This function aims to find weights \eqn{\hat\gamma_i(w)}, using logistic propensity scores,
#' that balance covariate means to a target vector, i.e.,
#' \deqn{\frac{1}{n} \sum_{i=1}^n \hat\gamma_i(W_i) X_i = \bar X_{\mathrm{target}}.}
#' With lasso regularization (`alpha = 1`), imbalance is controlled in the \eqn{\ell_\infty} sense,
#' allowing absolute slack of at most \eqn{\lambda} per covariate.
#'
#' For `target = "ATE"`, two logistic models are fit, one per arm, with
#' \deqn{\hat\gamma_i(1) = \frac{W_i}{\hat e^{(1)}(X_i)}, \quad
#'       \hat\gamma_i(0) = \frac{1 - W_i}{1 - \hat e^{(0)}(X_i)}, \quad
#'       \bar X_{\mathrm{target}} = \frac{1}{n} \sum_{i=1}^n X_i.}
#' \eqn{\hat e^{(w)}(X_i)} is the fitted propensity score for arm \eqn{w}.
#' For `target = "ATT"`, weights balance the control means:
#' \deqn{\hat\gamma_i = (1 - W_i) \frac{\hat e^{(0)}(X_i)}{1 - \hat e^{(0)}(X_i)}, \quad
#'       \bar X_{\mathrm{target}} = \frac{1}{\sum W_i} \sum_{i=1}^n W_i X_i.}
#'
#' @param X A numeric matrix or data frame with pre-treatment covariates.
#' @param W Treatment vector (0 = control, 1 = treated).
#' @param target The target estimand. Default is "ATE".
#' @param sample.weights Optional sample weights. If `NULL` (default), each unit receives the same weight.
#' @param max.imbalance Optional upper bound on the standardized covariate imbalance. For lasso penalization
#'   (`alpha = 1`), there is a one-to-one correspondence between the penalty parameter \eqn{\lambda} and
#'   the maximum allowable covariate imbalance. When supplied, `max.imbalance` is used to adjust the lambda
#'   sequence (via `lambda.min.ratio`) so that the generated sequence ends at the specified imbalance level.
#' @param nlambda Number of values for `lambda` if generated automatically. Default is 100.
#' @param lambda.min.ratio Ratio of smallest to largest lambda. Default is 1e-2.
#' @param lambda Optional `lambda` sequence. By default, it is constructed automatically using `nlambda`
#'   and `lambda.min.ratio` (or `max.imbalance`, if specified).
#' @param penalty.factor Penalty factor per feature. Default is 1 (i.e., each feature receives the same penalty).
#' @param groups Optional list of group indices for group penalization.
#' @param alpha Elastic net mixing parameter. Default is 1 (lasso), 0 corresponds to ridge.
#' @param standardize Whether to standardize the input matrix. Should only be `FALSE` if `X` already has
#'   zero-mean columns with unit variance. For `target = "ATT"`, standardization should be based on the treated group.
#' @param tol Coordinate descent convergence tolerance. Default is 1e-7.
#' @param maxit Maximum number of coordinate descent iterations. Default is 1e5.
#' @param verbose Whether to display information during fitting. Default is `FALSE`.
#' @param num.threads Number of threads to use. Default is 1.
#' @param ... Additional internal arguments passed to the solver.
#'
#' @return A fit balnet object.
#'
#' @references Sverdrup, Erik and Trevor Hastie.
#'  "balnet: Pathwise Estimation of Covariate Balancing Propensity Scores".
#'  arXiv preprint, arXiv:2602.18577, 2026.
#'
#' @examples
#' \donttest{
#' # Simulate data with confounding.
#' n <- 2000
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1.5 + exp(X[, 2] + X[, 3])))
#' Y <- W + 2 * log(1 + exp(X[, 1] + X[, 2] + X[, 3])) + rnorm(n)
#'
#' # Fit model targeting the ATE = E[Y(1)] - E[Y(0)].
#' # Two logistic models are fit: one for treated, one for control.
#' fit <- balnet(X, W, target = "ATE")
#'
#' # Print path summary.
#' print(fit)
#'
#' # Visualize the path.
#' plot(fit)
#'
#' # Plot the standardized covariate imbalance at given lambda.
#' # Note: lambda = 0 selects the final lambda in the sequence. Scalar values
#' # are applied to both arms.
#' plot(fit, lambda = 0)
#'
#' # Predict propensity scores at end of lambda path.
#' W.hat <- predict(fit, X, lambda = 0)
#'
#' # Get weights at end of lambda path.
#' ipw.weights <- ipw(fit, lambda = 0)
#'
#' # Estimate ATE using IPW weights.
#' mean(Y * (ipw.weights$treated - ipw.weights$control))
#' }
#'
#' @export
balnet <- function(
  X,
  W,
  target = c("ATE", "ATT", "treated", "control"),
  sample.weights = NULL,
  max.imbalance = NULL,
  nlambda = 100L,
  lambda.min.ratio = 1e-2,
  lambda = NULL,
  penalty.factor = NULL,
  groups = NULL,
  alpha = 1.0,
  standardize = TRUE,
  tol = 1e-7,
  maxit = as.integer(1e5),
  verbose = FALSE,
  num.threads = 1L,
  ...
)
{
  call <- match.call()
  target <- match.arg(target)
  if (is.matrix(X) || is.data.frame(X)) {
    if (!is.matrix(X) || !is.double(X)) {
      warning("X is not matrix of type 'double', an extra internal copy may occur.", immediate. = TRUE)
    }
    X <- as.matrix(X)
    if (!is.numeric(X) || anyNA(X) || any(0 %in% dim(X))) {
      stop("X should be numeric with no missing values.")
    }
  } else {
    stop("Invalid X input: should be a numeric matrix.")
  }
  if (!is.numeric(W) || length(W) != nrow(X) || anyNA(W) || any(W != 0 & W != 1)) {
    stop("W should be {0, 1} with length = nrow(X), with no missing values.")
  }
  if (!is.null(sample.weights) && (!is.numeric(sample.weights) || length(sample.weights) != nrow(X) || anyNA(sample.weights))) {
    stop("Invalid sample weights.")
  } else if (is.null(sample.weights)) {
    sample.weights <- rep_len(1, nrow(X))
  }
  if (is.character(standardize) && standardize == ".inplace") {
    inplace <- TRUE
    standardize <- TRUE
  } else if (is.logical(standardize)) {
    inplace <- FALSE
  } else {
    stop("Invalid standardize option.")
  }
  lambda <- validate_lambda(lambda)
  colnames <- if (is.null(colnames(X))) make.names(1:ncol(X)) else colnames(X)
  if (anyDuplicated(colnames)) {
    stop("Duplicate colnames detected in X.") # ensure colnames are valid as they are used later
  }
  validate_groups(groups, ncol(X), colnames)

  stan <- standardize(
    X,
    weights = if (target == "ATT") W * sample.weights else sample.weights,
    standardize = standardize,
    inplace = inplace,
    n_threads = num.threads
  )
  lambda.min.ratio <- get_lambda_min_ratio(lambda.min.ratio, max.imbalance, stan$X, W, sample.weights, target, alpha)

  if (target == "ATT") {
    target_scale = sum(sample.weights) / sum(sample.weights * W) # "n / n_1"
  } else {
    target_scale = 1
  }

  fit0 <- fit1 <- NULL
  lmdas0 <- lmdas1 <- NULL
  if (target %in% c("ATE", "ATT", "control")) {
    fit0 <- balnet.fit(
      stan = stan,
      y = 1 - W,
      weights = sample.weights,
      target_scale = target_scale,
      lambda = lambda[[1]],
      lmda_path_size = nlambda,
      min_ratio = lambda.min.ratio[[1]],
      penalty = penalty.factor,
      groups = groups,
      alpha = alpha,
      max_iters = maxit,
      tol = tol,
      progress_bar = verbose,
      progress_bar_prefix = if (target == "ATE") "Arm 0: " else "",
      n_threads = num.threads,
      ...
    )
    lmdas0 <- fit0$lmdas
  }
  if (target %in% c("ATE", "treated")) {
    fit1 <- balnet.fit(
      stan = stan,
      y = W,
      weights = sample.weights,
      target_scale = target_scale,
      lambda = lambda[[2]],
      lmda_path_size = nlambda,
      min_ratio = lambda.min.ratio[[2]],
      penalty = penalty.factor,
      groups = groups,
      alpha = alpha,
      max_iters = maxit,
      tol = tol,
      progress_bar = verbose,
      progress_bar_prefix = if (target == "ATE") "Arm 1: " else "",
      n_threads = num.threads,
      ...
    )
    lmdas1 <- fit1$lmdas
  }
  lambda <- list(control = lmdas0, treated = lmdas1)
  lambda.out <- lambda[!vapply(lambda, is.null, logical(1))]

  out <- list()
  class(out) <- "balnet"

  out[["call"]] <- call
  out[["X.dim"]] <- dim(X)
  out[["X.orig"]] <- X
  out[["W.orig"]] <- W
  out[["sample.weights"]] <- sample.weights
  out[["target"]] <- target
  out[["verbose"]] <- verbose
  out[["num.threads"]] <- num.threads
  out[["colnames"]] <- colnames
  out[["groups"]] <- groups
  out[["lambda"]] <- if (length(lambda.out) > 1) lambda.out else lambda.out[[1]]
  out[["_fit"]] <- list(control = fit0, treated = fit1)
  out[["_lambda"]] <- lambda

  out
}

#' Extract coefficients from a balnet object.
#'
#' @param object A `balnet` object.
#' @param lambda Value(s) of the penalty parameter `lambda` at which coefficients
#'   are required.
#'   * If `NULL` (default), the full lambda path from the fit is used.
#'   * If new values are supplied, linear interpolation is performed.
#'     For dual-arm fits (`target` = "ATE"), `lambda` can be a `list` or
#'     two-column `matrix`: the first element/column corresponds to the control
#'     arm and the second to the treatment.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Estimated logistic coefficients
#'  (for dual-arm fits, returns a list with entries for each arm).
#'
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATT model.
#' fit <- balnet(X, W, target = "ATT")
#'
#' # Extract coefficients.
#' coefs <- coef(fit)
#' }
#'
#' @method coef balnet
#' @export
coef.balnet <- function(
  object,
  lambda = NULL,
  ...
)
{
  lambda <- validate_lambda(lambda)

  coef0 <- coef1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    coef0 <- coef(object[["_fit"]]$control, lambda = lambda[[1]])
  }
  if (!is.null(object[["_fit"]]$treated)) {
    coef1 <- coef(object[["_fit"]]$treated, lambda = lambda[[2]])
  }
  out <- list(control = coef0, treated = coef1)
  out.nn <- out[!vapply(out, is.null, logical(1))]

  if (length(out.nn) > 1) {
    return(out.nn)
  } else {
    return(out.nn[[1]])
  }
}

#' Predict using a balnet object.
#'
#' @param object A `balnet` object.
#' @param newdata A numeric matrix.
#' @param lambda Value(s) of the penalty parameter `lambda` at which coefficients
#'   are required.
#'   * If `NULL` (default), the full lambda path from the fit is used.
#'   * If new values are supplied, linear interpolation is performed.
#'     For dual-arm fits (`target` = "ATE"), `lambda` can be a `list` or
#'     two-column `matrix`: the first element/column corresponds to the control
#'     arm and the second to the treatment.
#' @param type The type of predictions. Default is "response" (propensity scores).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Estimated predictions
#'  (for dual-arm fits, returns a list with entries for each arm).
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATT model.
#' fit <- balnet(X, W, target = "ATT")
#'
#' # Predict propensity scores.
#' W.hat <- predict(fit, X)
#' }
#'
#' @method predict balnet
#' @export
predict.balnet <- function(
  object,
  newdata,
  lambda = NULL,
  type = c("response"),
  ...
)
{
  lambda <- validate_lambda(lambda)
  type <- match.arg(type)
  if (missing(newdata)) {
    stop("newdata required for predictions.")
  }
  if (is.matrix(newdata) || is.data.frame(newdata)) {
    newdata <- as.matrix(newdata)
    if (!is.numeric(newdata) || anyNA(newdata) || ncol(newdata) != object[["X.dim"]][2]) {
      stop("newdata should be a numeric matrix with same columns as training data, with no missing values.")
    }
  } else {
    stop("Invalid newdata input: should be a numeric matrix.")
  }

  pred0 <- pred1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    pred0 <- 1 - predict(object[["_fit"]]$control, newdata, lambda = lambda[[1]], type = type)
  }
  if (!is.null(object[["_fit"]]$treated)) {
    pred1 <- predict(object[["_fit"]]$treated, newdata, lambda = lambda[[2]], type = type)
  }
  out <- list(control = pred0, treated = pred1)
  out.nn <- out[!vapply(out, is.null, logical(1))]

  dots <- list(...)
  .simplify <- is.null(dots[[".simplify"]])
  if (length(out.nn) > 1 || !.simplify) {
    return(out.nn)
  } else {
    return(out.nn[[1]])
  }
}

#' Print a balnet object.
#'
#' @param x A `balnet` object.
#' @param digits Number of digits to print.
#' @param max Total number of rows to show from the beginning and end of the path
#' @param ... Additional print arguments.
#'
#' @return Invisibly returns the printed information.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATT model.
#' fit <- balnet(X, W, target = "ATT")
#'
#' # Print path summary.
#' print(fit)
#' }
#'
#' @method print balnet
#' @export
print.balnet <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  max = 3,
  ...
)
{
  .print_compact <- function(txt) {
    if (length(txt) > 2 * max + 1) {
      out <- c(
        head(txt, max + 1),
        "...",
        tail(txt, max)
      )
    } else {
      out <- txt
    }
    cat(out, sep = "\n")
  }

  .get_metrics <- function(fit) {
    if (is.null(fit)) {
      out <- data.frame()
    } else {
      coeffs <- coef(fit)
      non.zero <- rowSums(coeffs$betas != 0)
      metric <- fit$devs
      lmdas <- fit$lmdas
      out <- data.frame(
        `Nonzero` = non.zero,
        `Mean |SMD|` = metric,
        Lambda = lmdas,
        check.names = FALSE
      )
    }
    out
  }

  cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n\n")
  df0 <- .get_metrics(x[["_fit"]]$control)
  df1 <- .get_metrics(x[["_fit"]]$treated)
  to.print <- rbind(df0, df1)
  attr(to.print, "row.names") <- c(rownames(df0), rownames(df1))
  txt <- utils::capture.output(print(to.print, digits = digits, ...)) # align the printout if two arms

  if (!is.null(x[["_fit"]]$control)) {
    fit <- x[["_fit"]]$control
    cat("Control ", "(path: ", nrow(df0), "/", length(fit$lmda_path), ")\n", sep = "")
    .print_compact(txt[1:(nrow(df0) + 1)])
    txt <- txt[-(2:(nrow(df0) + 1))]
    if (!is.null(x[["_fit"]]$treated)) cat("\n")
  }
  if (!is.null(x[["_fit"]]$treated)) {
    fit <- x[["_fit"]]$treated
    cat("Treated ", "(path: ", nrow(df1), "/", length(fit$lmda_path), ")\n", sep = "")
    .print_compact(txt)
  }
  out <- list(control = df0, treated = df1)
  out.nn <- out[vapply(out, length, integer(1)) > 0]

  dots <- list(...)
  .simplify <- is.null(dots[[".simplify"]])
  if (length(out.nn) > 1 || !.simplify) {
    invisible(out.nn)
  } else {
    invisible(out.nn[[1]])
  }
}

#' Plot diagnostics for a `balnet` object.
#'
#' Shows effective sample size (ESS) and percent bias reduction (PBR; reduction
#' in mean absolute imbalance) along the regularization path, computed from IPW
#' weights and normalized to percentages. The right-hand axis maps these values
#' to the coefficient of variation (CV) of the weights.
#' Supplying the `lambda` argument displays the standardized covariate imbalance
#' \eqn{(\bar X_{\mathrm{weighted}} - \bar X_{\mathrm{target}}) / \sigma_{\mathrm{target}}},
#' computed using the IPW weights at the specified `lambda`.
#'
#' @param x A `balnet` object.
#' @param lambda If NULL (default) diagnostics over the lambda path is shown.
#' Otherwise, covariate balance at provided lambda value is shown
#' (if target = "ATE", lambda can be a 2-vector, arm 0 and arm 1.)
#' @param groups Optional named list of contiguous covariate index ranges to
#'   aggregate into a single variable before computing covariate imbalance
#'   (e.g., `list(demographics = 4:12)`).
#' @param max The number of covariates to display in covariate balance plot. Defaults to all covariates.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the information underlying the plot.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATT model.
#' fit <- balnet(X, W, target = "ATT")
#'
#' # Plot the five covariates with the largest unweighted imbalance
#' plot(fit, lambda = 0, max = 5)
#' }
#'
#' @method plot balnet
#' @export
plot.balnet <- function(
  x,
  lambda = NULL,
  groups = NULL,
  max = NULL,
  ...
)
{
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))

  if (is.null(lambda)) {
    plot_func <- plot_path
    get_metrics <- get_path
  } else {
    if (length(lambda) > 2) {
      stop("Can only plot for a single value of lambda.")
    }
    plot_func <- plot_smd
    get_metrics <- get_smd
    if (length(lambda) == 1) {
      lambda <- list(c(Inf, lambda[[1]]), c(Inf, lambda[[1]]))
    } else {
      lambda <- list(c(Inf, lambda[[1]]), c(Inf, lambda[[2]]))
    }
  }

  lambda.orig <- x[["_lambda"]]
  W.orig <- x[["W.orig"]]
  W.hat <- predict.balnet(x, x[["X.orig"]], lambda = lambda, type = "response", .simplify = FALSE)

  stats0 <- stats1 <- NULL
  if (!is.null(x[["_fit"]]$control)) {
    stats0 <- get_metrics(
      x,
      1 - W.hat$control,
      1 - W.orig,
      lambda = lambda.orig$control,
      groups = groups,
      devs = x[["_fit"]]$control$devs
    )
    if (!is.null(x[["_fit"]]$treated)) graphics::par(mfrow = c(1, 2))
    plot_func(stats0, max = max)
    if (x[["target"]] == "ATE") graphics::mtext("Control", side = 3, line = 1, adj = 0)
  }

  if (!is.null(x[["_fit"]]$treated)) {
    stats1 <- get_metrics(
      x,
      W.hat$treated,
      W.orig,
      lambda = lambda.orig$treated,
      devs = x[["_fit"]]$treated$devs,
      groups = groups
    )
    plot_func(stats1, max = max)
    if (x[["target"]] == "ATE") graphics::mtext("Treated", side = 3, line = 1, adj = 0)
  }
  out <- list(control = stats0, treated = stats1)
  out.nn <- out[!vapply(out, is.null, logical(1))]

  if (length(out.nn) > 1) {
    invisible(out.nn)
  } else {
    invisible(out.nn[[1]])
  }
}

#' Extract IPW weights from a balnet object.
#'
#' @param object A `balnet` object.
#' @param lambda Value(s) of the penalty parameter `lambda` at which coefficients
#'   are required.
#'   * If `NULL` (default), the full lambda path from the fit is used.
#'   * If new values are supplied, linear interpolation is performed.
#'     For dual-arm fits (`target` = "ATE"), `lambda` can be a `list` or
#'     two-column `matrix`: the first element/column corresponds to the control
#'     arm and the second to the treatment.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Estimated IPW weights.
#'  (for contrast fits, `target` = "ATE" or "ATT", returns a list with entries for each arm).
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATT model.
#' fit <- balnet(X, W, target = "ATT")
#'
#' # Extract IPW weights.
#' wts <- ipw(fit, lambda = 0)
#' }
#'
#' @export
ipw <- function(
  object,
  lambda = NULL,
  ...
)
{
  UseMethod("ipw")
}

#' @rdname ipw
#' @method ipw balnet
#' @export
ipw.balnet <- function(
  object,
  lambda = NULL,
  ...
)
{
  lambda <- validate_lambda(lambda)
  X.orig <- object[["X.orig"]]
  W.orig <- object[["W.orig"]]
  sample.weights <- object[["sample.weights"]]
  target <- object[["target"]]

  W.hat <- predict.balnet(object, X.orig, lambda = lambda, type = "response", .simplify = FALSE)

  ipw0 <- ipw1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    ipw0 <- matrix(0, nrow = nrow(W.hat$control), ncol = ncol(W.hat$control))
    ipw0[W.orig == 0, ] <- 1 / (1 - W.hat$control[W.orig == 0, ])
    if (target == "ATT") {
      ipw0[W.orig == 0, ] <- ipw0[W.orig == 0, ] * W.hat$control[W.orig == 0, ]
      ipw1 <- matrix(0, nrow = nrow(W.hat$control), ncol = ncol(W.hat$control))
      ipw1[W.orig == 1, ] <- sample.weights[W.orig == 1]
    }
    ipw0 <- ipw0 * sample.weights
  }
  if (!is.null(object[["_fit"]]$treated)) {
    ipw1 <- matrix(0, nrow = nrow(W.hat$treated), ncol = ncol(W.hat$treated))
    ipw1[W.orig == 1, ] <- 1 / W.hat$treated[W.orig == 1, ]
    ipw1 <- ipw1 * sample.weights
  }
  out <- list(control = ipw0, treated = ipw1)
  out.nn <- out[!vapply(out, is.null, logical(1))]

  if (length(out.nn) > 1) {
    return(out.nn)
  } else {
    return(out.nn[[1]])
  }
}

get_path <- function(fit, W.hat, W, ..., lambda, devs) {
  target <- fit[["target"]]
  sample.weights <- fit[["sample.weights"]]

  ipw <- matrix(0, nrow = nrow(W.hat), ncol = ncol(W.hat))
  ipw[W == 1, ] <- 1 / W.hat[W == 1, ]
  if (target == "ATT") {
    ipw <- (1 - W.hat) * ipw
  }
  ess <- (colSums(sample.weights * ipw)^2 / colSums(sample.weights * ipw^2)) /
    sum(W * sample.weights) * 100
  pbr <- (1 - devs / devs[1]) * 100 # "devs" stores mean abs(SMD)

  data.frame(lambda = lambda, ess = ess, pbr = pbr)
}

get_smd <- function(fit, W.hat, W, ..., groups = NULL) {
  target <- fit[["target"]]
  sample.weights <- fit[["sample.weights"]]
  X <- fit[["X.orig"]]
  colnames <- fit[["colnames"]]
  # if groups present, we calculate SMDs on group-level means
  if (!is.null(groups)) {
    X <- collapse_X(X, groups, colnames)
    colnames <- colnames(X)
  }

  # ATT SMD: (\weighted \bar X_C - \bar X_T) / S_T
  if (target == "ATT") {
    X.stats <- col_stats(X, fit[["W.orig"]] * sample.weights, compute_sd = TRUE)
  } else {
    X.stats <- col_stats(X, sample.weights, compute_sd = TRUE)
  }
  X.stats$scale[X.stats$scale <= 0] <- 1

  ipw <- matrix(0, nrow = nrow(W.hat), ncol = ncol(W.hat))
  ipw[W == 1, ] <- 1 / W.hat[W == 1, ]
  if (target == "ATT") {
    ipw <- (1 - W.hat) * ipw
  }

  smd <- col_stats(X, ipw * sample.weights, n_threads = fit[["num.threads"]])$center
  smd <- sweep(smd, 2L, X.stats$center, `-`, check.margin	= FALSE)
  smd <- sweep(smd, 2L, X.stats$scale, `/`, check.margin = FALSE)

  out <- data.frame(t(smd), row.names = colnames)
  colnames(out)[1:2] <- c("lambda.max", "lambda")

  out
}

plot_path <- function(stats, ...) {
  lambdas <- stats$lambda
  pbr <- stats$pbr
  ess <- stats$ess

  graphics::plot(lambdas[lambdas > 0], pbr[lambdas > 0],
    log = "x",
    type = "l",
    xlim = rev(range(lambdas[lambdas > 0])),
    ylim = c(min(0, min(pbr)), 100),
    xlab = expression(Log(lambda)),
    ylab = "Percent  (CV of weights)"
  )
  graphics::points(lambdas[lambdas > 0], ess[lambdas > 0], type = "l", col = "dodgerblue3")
  graphics::mtext("PBR", side = 3, adj = 1, line = 1 )
  graphics::mtext("ESS", side = 3, adj = 1, line = 0, col = "dodgerblue3")
  graphics::abline(h = 0)
  # coeff of variation of weights is sqrt(100 / ESS%  - 1)
  cv.ticks <- c(0, 0.5, 1, 1.5, 2, 3)
  ess.pos <- 100 / (1 + cv.ticks^2)
  graphics::axis(side = 4, at = ess.pos, labels = cv.ticks, las = 0, cex.axis = 0.8)
}

plot_smd <- function(stats, ..., max = NULL) {
    labels <- rownames(stats)
    if (is.null(max)) {
      max <- length(labels)
    }
    max <- min(max, length(labels))
    order <- order(abs(stats$lambda.max), decreasing = TRUE)
    display.idx <- rev(order[1:max])

    graphics::plot(
      stats$lambda.max[display.idx],
      1:max,
      xlim = c(min(-0.1, min(stats$lambda.max)), max(0.1, max(stats$lambda.max))),
      xlab = "Standardized mean diff.",
      ylab = "",
      pch = 20,
      yaxt = "n"
    )
    graphics::axis(2, at = 1:max, labels = labels[display.idx], las = 1, cex.axis = 0.7)
    graphics::points(stats$lambda[display.idx], 1:max, pch = 20, col = "dodgerblue3")
    graphics::abline(v = 0)
    graphics::abline(v = c(-0.1, 0.1), lty = 2, col = "gray70")
    graphics::mtext(expression(lambda^{max}), side = 3, adj = 1, line = 1 )
    graphics::mtext(expression(lambda^{phantom(max)}), side = 3, adj = 1, line = 0, col = "dodgerblue3")
}
