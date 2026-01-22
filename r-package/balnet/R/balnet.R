#' Pathwise estimation of covariate balancing propensity scores.
#'
#' @param X Feature matrix.
#' @param W Treatment vector (0: control, 1: treated).
#' @param target The target estimand. Default is ATE.
#' @param sample.weights Optional sample weights. If `NULL` (default), then each unit receives the same weight.
#' @param nlambda Number of values for `lambda`, if generated automatically. Default is 100.
#' @param lambda.min.ratio Ratio between smallest and largest value of lambda. Default is 1e-2.
#' @param lambda Optional `lambda` sequence. By default, the `lambda` sequence is constructed automatically using `nlambda` and `lambda.min.ratio`.
#' @param penalty.factor Penalty factor per feature. Default is 1 (i.e, each feature recieves the same penalty).
#' @param groups A list of group indices.
#' @param alpha Elastic net mixing parameter. Default is 1 (lasso). 0 is ridge.
#' @param standardize Whether to standardize the input matrix. This should only be set to `FALSE` if `X` already has zero-mean columns with unit variances
#'  (for `target = "ATT"`, `X` is expected to be standardized with the treated information).
#'  It can also be set to "inplace", which overwrites the input matrix `X` with standardized information.
#' @param thresh Coordinate descent convergence tolerance, default 1e-7.
#' @param maxit Maximum total number of coordinate descent iterations, default is 1e5.
#' @param verbose Whether to display information during fitting. Default is `FALSE`.
#' @param num.threads Number of threads, default is 1.
#' @param ... Additional internal arguments passed to solver.
#'
#' @return A fit balnet object.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' fit <- balnet(X, W)
#'
#' # Print path summary.
#' print(fit)
#'
#' # Plot path diagnostics.
#' plot(fit)
#'
#' # Plot covariate imbalance at the end of the path (closest to lambda = 0).
#' plot(fit, lambda = 0)
#'
#' # Predict propensity scores.
#' pp <- predict(fit, X)
#'
#' # Extract coefficients.
#' coefs <- coef(fit)
#' }
#'
#' @export
balnet <- function(
  X,
  W,
  target = c("ATE", "ATT", "treated", "control"),
  sample.weights = NULL,
  nlambda = 100L,
  lambda.min.ratio = 1e-2,
  lambda = NULL,
  penalty.factor = NULL,
  groups = NULL,
  alpha = 1.0,
  standardize = TRUE,
  thresh = 1e-7,
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
    if (!is.numeric(X) || anyNA(X)) {
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
  lambda.in <- validate_lambda(lambda)
  if (is.character(standardize) && standardize == "inplace") {
    inplace <- TRUE
    standardize <- TRUE
  } else if (is.logical(standardize)) {
    inplace <- FALSE
  } else {
    stop("Invalid standardize option.")
  }
  colnames <- if (is.null(colnames(X))) make.names(1:ncol(X)) else colnames(X)
  validate_groups(groups, ncol(X), colnames)

  stan <- standardize(
    X,
    weights = if (target == "ATT") W * sample.weights else sample.weights,
    standardize = standardize,
    inplace = inplace,
    n_threads = num.threads
  )
  if (target == "ATT") {
    target_scale = sum(sample.weights) / sum(sample.weights * W) # "n / n_1"
  } else {
    target_scale = 1
  }

  fit0 <- fit1 <- NULL
  lmdas0 <- lmdas1 <- NULL
  if (target %in% c("ATE", "ATT", "control")) {
    if (verbose) message("Fitting arm: 0")
    fit0 <- balnet.fit(
      stan = stan,
      y = 1 - W,
      weights = sample.weights,
      target_scale = target_scale,
      lambda = lambda.in[[1]],
      lmda_path_size = nlambda,
      min_ratio = lambda.min.ratio,
      penalty = penalty.factor,
      groups = groups,
      alpha = alpha,
      max_iters = maxit,
      tol = thresh,
      progress_bar = verbose,
      n_threads = num.threads,
      ...
    )
    lmdas0 <- fit0$lmdas
  }
  if (target %in% c("ATE", "treated")) {
    if (verbose) message("Fitting arm: 1")
    fit1 <- balnet.fit(
      stan = stan,
      y = W,
      weights = sample.weights,
      target_scale = target_scale,
      lambda = lambda.in[[2]],
      lmda_path_size = nlambda,
      min_ratio = lambda.min.ratio,
      penalty = penalty.factor,
      groups = groups,
      alpha = alpha,
      max_iters = maxit,
      tol = thresh,
      progress_bar = verbose,
      n_threads = num.threads,
      ...
    )
    lmdas1 <- fit1$lmdas
  }
  lambda <- list(control = lmdas0, treated = lmdas1)

  out <- list()
  class(out) <- "balnet"

  out[["call"]] <- call
  out[["X.orig"]] <- X
  out[["W.orig"]] <- W
  out[["sample.weights"]] <- sample.weights
  out[["target"]] <- target
  out[["verbose"]] <- verbose
  out[["num.threads"]] <- num.threads
  out[["colnames"]] <- colnames
  out[["groups"]] <- groups
  out[["lambda"]] <- lambda[sapply(lambda, length) > 0]
  out[["_fit"]] <- list(control = fit0, treated = fit1)
  out[["_cache"]] <- new.env(parent = emptyenv()) # saves plotting statistics etc.

  out
}

#' Extract coefficients from a balnet object.
#'
#' @param object A `balnet` object.
#' @param lambda Value(s) of the penalty parameter `lambda` at which coefficients
#'   are required.
#'   * If `NULL` (default), the full lambda path from the fit is used
#'    (if new values are supplied, linear interpolation is performed).
#'   * For dual-arm fits (control and treatment), `lambda` can be a `list` or
#'     two-column `matrix`: the first element/column corresponds to the control
#'     arm and the second to the treatment arm.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The estimated coefficients.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' fit <- balnet(X, W)
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
  lambda.in <- validate_lambda(lambda)

  coef1 <- coef0 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    coef0 <- coef(object[["_fit"]]$control, lambda = lambda.in[[1]])
  }
  if (!is.null(object[["_fit"]]$treated)) {
    coef1 <- coef(object[["_fit"]]$treated, lambda = lambda.in[[2]])
  }
  out <- list(control = coef0, treated = coef1)

  out[sapply(out, length) > 0]
}

#' Predict using a balnet object.
#'
#' @param object A `balnet` object.
#' @param newx A numeric matrix.
#' @param lambda Value(s) of the penalty parameter `lambda` at which coefficients
#'   are required.
#'   * If `NULL` (default), the full lambda path from the fit is used
#'    (if new values are supplied, linear interpolation is performed).
#'   * For dual-arm fits (control and treatment), `lambda` can be a `list` or
#'     two-column `matrix`: the first element/column corresponds to the control
#'     arm and the second to the treatment arm.
#' @param type The type of predictions.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Predictions.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' fit <- balnet(X, W)
#'
#' # Predict propensity scores.
#' pp <- predict(fit, X)
#' }
#'
#' @method predict balnet
#' @export
predict.balnet <- function(
  object,
  newx,
  lambda = NULL,
  type = c("response"),
  ...
)
{
  lambda.in <- validate_lambda(lambda)
  type <- match.arg(type)
  dot.args <- list(...)
  if (missing(newx) && !is.null(dot.args$newdata)) {
    newx <- dot.args$newdata
  } else if (missing(newx)) {
    stop("newx required for predictions.")
  }
  if ("s" %in% names(dot.args) && is.null(lambda)) {
    lambda <- dot.args$s
  }
  if (is.matrix(newx) || is.data.frame(newx)) {
    newx <- as.matrix(newx)
    if (!is.numeric(newx) || anyNA(newx) || ncol(newx) != ncol(object[["X.orig"]])) {
      stop("X should be a numeric with same columns as training data, with no missing values.")
    }
  } else {
    stop("Invalid X input: should be a numeric matrix.")
  }

  pred0 <- pred1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    pred0 <- predict(object[["_fit"]]$control, newx, lambda = lambda.in[[1]], type = type)
  }
  if (!is.null(object[["_fit"]]$treated)) {
    pred1 <- predict(object[["_fit"]]$treated, newx, lambda = lambda.in[[2]], type = type)
  }
  out <- list(control = pred0, treated = pred1)

  out[sapply(out, length) > 0]
}

#' Print a balnet object.
#'
#' @param x A `balnet` object.
#' @param digits Number of digits to print.
#' @param max Total number of rows to show from the beginning and end of the path
#' @param ... Additional print arguments.
#'
#' @return Invisibly returns a data.frame with the printed information.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' fit <- balnet(X, W)
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

  invisible(out[sapply(out, length) > 0])
}

#' Plot diagnostics for a `balnet` object.
#'
#' @param x A `balnet` object.
#' @param lambda If NULL (default) diagnostics over the lambda path is shown.
#' Otherwise, diagnostics for a single lambda value is shown.
#' (if target = "ATE", lambda can be a 2-vector, arm 0 and arm 1.)
#' @param groups A list of group indices.
#' @param max The number of covariates to display in balance plot. Defaults to all covariates.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns a list with the information underlying the plot.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' fit <- balnet(X, W)
#'
#' # Plot path diagnostics.
#' plot(fit)
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

  if (!is.null(lambda) && length(lambda) > 2) {
    stop("Can only plot for a single value of lambda.")
  }
  if (length(lambda) == 1L) {
    lambda.in <- list(lambda[1], lambda[1])
  } else {
    lambda.in <- lambda
  }
  plot_func <- if (is.null(lambda)) `plot_path` else `plot_smd`

  if (!is.null(groups)) {
    rm(list = intersect(c("stats0", "stats1"), ls(envir = x[["_cache"]])), envir = x[["_cache"]])
  }

  lambdas <- x[["lambda"]]
  W.orig <- x[["W.orig"]]
  pp <- get0("pp", envir = x[["_cache"]], inherits = FALSE)

  if (is.null(pp)) {
    pp <- predict(x, x[["X.orig"]], lambda = NULL, type = "response")
    x[["_cache"]]$pp <- pp
  }

  stats0 <- get0("stats0", envir = x[["_cache"]], inherits = FALSE)
  stats1 <- get0("stats1", envir = x[["_cache"]], inherits = FALSE)
  if (!is.null(x[["_fit"]]$control)) {
    if (is.null(stats0)) {
      stats0 <- get_metrics(lambdas$control, pp$control, 1 - W.orig, groups, x)
      if (is.null(groups)) x[["_cache"]]$stats0 <- stats0
    }
    if (!is.null(x[["_fit"]]$treated)) {
      graphics::par(mfrow = c(1, 2))
    }
    plot_func(stats0, lambda.in[[1]], max)
    if (x[["target"]] == "ATE") graphics::mtext("Control", side = 3, line = 1, adj = 0)
  }

  if (!is.null(x[["_fit"]]$treated)) {
    if (is.null(stats1)) {
      stats1 <- get_metrics(lambdas$treated, pp$treated, W.orig, groups, x)
      if (is.null(groups)) x[["_cache"]]$stats1 <- stats1
    }
    plot_func(stats1, lambda.in[[2]], max)
    if (x[["target"]] == "ATE") graphics::mtext("Treated", side = 3, line = 1, adj = 0)
  }
  out <- list(control = stats0, treated = stats1)

  invisible(out[sapply(out, length) > 0])
}

get_metrics <- function(lambdas, pp, W, groups, fit) {
  target <- fit[["target"]]
  X <- fit[["X.orig"]]
  colnames <- fit[["colnames"]]
  # if groups present, we calculate SMDs on group-level means
  if (!is.null(groups)) {
    X <- collapse_X(X, groups, colnames)
    colnames <- colnames(X)
  }

  # ATT SMD: (\weighted \bar X_C - \bar X_T) / S_T
  if (target == "ATT") {
    X.stats <- col_stats(X, fit[["W.orig"]], compute_sd = TRUE)
  } else {
    X.stats <- col_stats(X, compute_sd = TRUE)
  }
  X.stats$scale[X.stats$scale <= 0] <- 1

  ipw <- W  / pp
  if (target == "ATT") {
    ipw <- (1 - pp) * ipw
  }
  ess <- (colSums(ipw)^2 / colSums(ipw^2)) / sum(W) * 100

  smd <- col_stats(X, ipw, n_threads = fit[["num.threads"]])$center
  smd <- sweep(smd, 2L, X.stats$center, `-`, check.margin	= FALSE)
  smd <- sweep(smd, 2L, X.stats$scale, `/`, check.margin = FALSE)
  pbr <- (1 - rowSums(abs(smd)) / sum(abs(smd[1, ]))) * 100
  colnames(smd) <- colnames

  list(
    pth = cbind(lambda = lambdas, ess = ess, pbr = pbr),
    smd = cbind(lambda = lambdas, smd)
  )
}

plot_path <- function(stats, lambda, ...) {
  lambdas <- stats[["pth"]][, "lambda"]
  pbr <- stats[["pth"]][, "pbr"]
  ess <- stats[["pth"]][, "ess"]

  graphics::plot(lambdas[lambdas > 0], pbr[lambdas > 0],
    log = "x",
    type = "l",
    xlim = rev(range(lambdas[lambdas > 0])),
    ylim = c(min(0, min(pbr)), 100),
    xlab = expression(Log(lambda)),
    ylab = "percent"
  )
  graphics::points(lambdas[lambdas > 0], ess[lambdas > 0], type = "l", col = "dodgerblue3")
  graphics::mtext("PBR", side = 3, adj = 1, line = 1 )
  graphics::mtext("ESS", side = 3, adj = 1, line = 0, col = "dodgerblue3")
  graphics::abline(h = 0)
}

plot_smd <- function(stats, lambda, max = NULL, ...) {
    lambdas <- stats[["smd"]][, "lambda"]
    smd <- stats[["smd"]][, -1]
    labels <- colnames(smd)
    if (is.null(max)) {
      max <- length(labels)
    }
    max <- min(max, length(labels))
    order <- order(abs(smd[1, ]), decreasing = TRUE)
    display.idx <- rev(order[1:max])
    lmda.ix <- max(findInterval(-lambda, -lambdas), 1) # find nearest lambda on decreasing path

    graphics::plot(
      smd[1, display.idx],
      1:max,
      xlim = c(min(-0.1, min(smd[1, ])), max(0.1, max(smd[1, ]))),
      xlab = "Standardized mean diff.",
      ylab = "",
      pch = 20,
      yaxt = "n"
    )
    graphics::axis(2, at = 1:max, labels = labels[display.idx], las = 1, cex.axis = 0.7)
    graphics::points(smd[lmda.ix, display.idx], 1:max, pch = 20, col = "dodgerblue3")
    graphics::abline(v = 0)
    graphics::abline(v = c(-0.1, 0.1), lty = 2, col = "gray70")
    graphics::mtext(expression(lambda^{max}), side = 3, adj = 1, line = 1 )
    graphics::mtext(expression(lambda^{phantom(max)}), side = 3, adj = 1, line = 0, col = "dodgerblue3")
}
