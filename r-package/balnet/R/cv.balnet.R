#' Tuning via cross-validation for balnet.
#'
#' @param X A numeric matrix or data frame with pre-treatment covariates.
#' @param W Treatment vector (0: control, 1: treated).
#' @param type.measure The loss to minimize for cross-validation.
#'  Default is balance loss (e.g., Zhiqiang (2020)).
#'  For "imbalance", the criterion is mean covariate imbalance (e.g., Zhao (2019)).
#' @param nfolds The number of folds used for cross-validation, default is 10.
#' @param foldid An optional `n`-vector specifying which fold 1 to `nfold` a sample belongs to.
#' If NULL, this defaults to `sample(rep(seq(nfolds), length.out = nrow(X)))`.
#' @param ... Arguments for \code{\link{balnet}}.
#'
#' @return A fit cv.balnet object.
#'
#' @references Tan, Zhiqiang.
#'  "Regularized calibrated estimation of propensity scores with model misspecification and high-dimensional data."
#'  Biometrika 107(1), 2020.
#' @references Zhao, Qingyuan.
#'  "Covariate balancing propensity score by tailored loss functions."
#'  Annals of Statistics 47 (2), 2019.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' cv.fit <- cv.balnet(X, W)
#'
#' # Print CV summary.
#' print(cv.fit)
#'
#' # Plot at cross-validated lambda.
#' plot(cv.fit)
#'
#' # Predict at cross-validated lambda.
#' W.hat <- predict(cv.fit, X)
#'
#' }
#'
#' @export
cv.balnet <- function(
  X,
  W,
  type.measure = c("balance.loss", "imbalance"),
  nfolds = 10,
  foldid = NULL,
  ...
)
{
  dot.args <- list(...)
  type.measure <- match.arg(type.measure)
  X.stats <- NULL
  if (type.measure == "balance.loss") {
    get_loss <- get_balance_loss
  } else if (type.measure == "imbalance") {
    get_loss <- get_imbalance
    X.stats <- col_stats(X, dot.args[["sample.weights"]], compute_sd = TRUE)
    X.stats$scale[X.stats$scale <= 0] <- 1
  }

  if (is.null(foldid)) {
    nfolds <- max(nfolds, 3)
    foldid <- sample(rep(seq(nfolds), length.out = nrow(X)))
  } else {
    if (length(foldid) != length(W)) {
      stop("Invalid `foldid`.")
    }
    nfolds <- max(foldid)
  }

  if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) cat("Fitting full model\n")
  fit.full <- balnet(X, W, ...)
  lambda.full <- fit.full[["_lambda"]]
  sample.weights <- fit.full[["sample.weights"]]

  cv.list <- list()
  for (k in 1:nfolds) {
    if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) cat(sprintf("\nFold: %d/%d\n", k, nfolds))
    test <- foldid == k
    train <- !test
    X.train <- X[train, , drop = FALSE]
    W.train <- W[train]
    dot.args[["sample.weights"]] <- sample.weights[train]
    fit.train <- do.call(balnet, c(list(X = X.train, W = W.train, standardize = ".inplace"), dot.args))

    X.test <- X[test, , drop = FALSE]
    W.test <- W[test]
    sample.weights.test <- sample.weights[test]
    loss <- do.call(get_loss, list(fit.train, X.test, W.test, sample.weights.test, lambda.full, X.stats = X.stats))
    cv.list[[k]] <- loss
  }
  cv.mean0 <- cv.mean1 <- NULL
  idx.min0 <- idx.min1 <- NULL
  lambda.min0 <- lambda.min1 <- NULL
  if (!is.null(cv.list[[1]][["control"]])) {
    cv.mean0 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, "control")),
                                nrow = length(cv.list),
                                ncol = length(lambda.full$control),
                                byrow = TRUE))
    idx.min0 <- which.min(cv.mean0)
    lambda.min0 <- lambda.full[["control"]][idx.min0]
  }
  if (!is.null(cv.list[[1]][["treated"]])) {
    cv.mean1 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, "treated")),
                                nrow = length(cv.list),
                                ncol = length(lambda.full$treated),
                                byrow = TRUE))
    idx.min1 <- which.min(cv.mean1)
    lambda.min1 <- lambda.full[["treated"]][idx.min1]
  }
  lambda.min <- list(control = lambda.min0, treated = lambda.min1)
  lambda.min.out <- lambda.min[!vapply(lambda.min, is.null, logical(1))]
  cv.info <- list(
    "cv.mean" = list(control = cv.mean0, treated = cv.mean1),
    "idx.min" = list(control = idx.min0, treated = idx.min1),
    "lambda.min" = lambda.min,
    "type.measure" = type.measure
  )

  fit.full[["lambda.min"]] <- if (length(lambda.min.out) > 1) lambda.min.out else lambda.min.out[[1]]
  fit.full[["_cv.info"]] <- cv.info
  fit.full[["call"]] <- match.call()
  class(fit.full) <- c("cv.balnet", class(fit.full))

  fit.full
}

#' Extract coefficients from a cv.balnet object.
#'
#' @param object A `cv.balnet` object.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Estimated logistic coefficients
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
#' cv.fit <- cv.balnet(X, W, target = "ATT")
#'
#' # Extract coefficients at cross-validated lambda.
#' coefs <- coef(cv.fit)
#' }
#'
#' @method coef cv.balnet
#' @export
coef.cv.balnet <- function(
  object,
  lambda = "lambda.min",
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- object[["_cv.info"]]$lambda.min
  }

  coef.balnet(object, lambda = lambda)
}

#' Predict using a cv.balnet object.
#'
#' @param object A `cv.balnet` object.
#' @param newdata A numeric matrix.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
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
#' cv.fit <- cv.balnet(X, W, target = "ATT")
#'
#' # Predict propensity scores at cross-validated lambda.
#' W.hat <- predict(cv.fit, X)
#' }
#'
#' @method predict cv.balnet
#' @export
predict.cv.balnet <- function(
  object,
  newdata,
  lambda = "lambda.min",
  type = c("response"),
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- object[["_cv.info"]]$lambda.min
  }

  predict.balnet(object, newdata, lambda = lambda, type = type)
}

#' Print a cv.balnet object.
#'
#' @param x A `cv.balnet` object.
#' @param digits Number of digits to print.
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
#' cv.fit <- cv.balnet(X, W, target = "ATT")
#'
#' # Print CV summary.
#' print(cv.fit)
#' }
#'
#' @method print cv.balnet
#' @export
print.cv.balnet <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  ...
)
{
  cat("Call: ", paste(deparse(x$call), collapse = "\n"), "\n\n")

  utils::capture.output(out <- print.balnet(x, digits = digits, .simplify = FALSE, ...))
  df0 <- df1 <- data.frame()
  if (!is.null(x[["_fit"]]$control)) {
    idx.min0 <- x[["_cv.info"]]$idx.min[[1]]
    df0 <- cbind(Arm = "Control", out$control[idx.min0, ], Index = idx.min0)
  }
  if (!is.null(x[["_fit"]]$treated)) {
    idx.min1 <- x[["_cv.info"]]$idx.min[[2]]
    df1 <- cbind(Arm = "Treated", out$treated[idx.min1, ], Index = idx.min1)
  }

  type.measure <- paste0("type.measure = ", x[["_cv.info"]]$type.measure)
  cat("Cross-validated lambda minimizing ", type.measure, ":\n", sep = "")
  print(rbind(df0, df1), digits = digits, row.names = FALSE, right = FALSE)
}

#' Summarize a cv.balnet object.
#'
#' @param object `cv.balnet` object.
#' @param ... Additional summary arguments.
#'
#' @return Returns the printed information.
#'
#' @method summary cv.balnet
#' @export
summary.cv.balnet <- function(object, ...) {
  print(object, ...)
}

#' Plot diagnostics for a `cv.balnet` object.
#'
#' @param x A `cv.balnet` object.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
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
#' cv.fit <- cv.balnet(X, W, target = "ATT")
#'
#' # Plot at cross-validated lambda.
#' plot(cv.fit)
#' }
#'
#' @method plot cv.balnet
#' @export
plot.cv.balnet <- function(
  x,
  lambda = "lambda.min",
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- x[["_cv.info"]]$lambda.min
  }

  plot.balnet(x, lambda = lambda, ...)
}

#' @rdname balweights
#' @method balweights cv.balnet
#' @export
balweights.cv.balnet <- function(
  object,
  lambda = "lambda.min",
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- object[["_cv.info"]]$lambda.min
  }

  balweights.balnet(object, lambda = lambda)
}

get_balance_loss <- function(object, X.test, W.test, sample.weights, lambda, ...) {
  .balance_loss <- function(W, eta) {
    colSums(sample.weights * (W * exp(-eta) + (1 - W) * eta)) / sum(sample.weights)
  }

  lambda <- validate_lambda(lambda)
  eta <- predict(object, X.test, lambda = lambda, type = "link", .simplify = FALSE)

  loss0 <- loss1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    loss0 <- .balance_loss(1 - W.test, eta$control)
  }
  if (!is.null(object[["_fit"]]$treated)) {
    loss1 <- .balance_loss(W.test, eta$treated)
  }
  out <- list(control = loss0, treated = loss1)

  out[!vapply(out, is.null, logical(1))]
}

get_imbalance <- function(object, X.test, W.test, sample.weights, lambda, X.stats, ...) {
  .imbalance <- function(W, W.hat) {
    # held-out balancing p-scores might be zero, so we need to clip.
    # (not as big a concern in the held out balance loss)
    W.hat[W == 1, ] <- pmax(W.hat[W == 1, ], 0.001)

    ipw <- matrix(0, nrow = nrow(W.hat), ncol = ncol(W.hat))
    ipw[W == 1, ] <- 1 / W.hat[W == 1, ]

    smd <- col_stats(X.test, ipw * sample.weights, n_threads = object[["num.threads"]])$center
    smd <- sweep(smd, 2L, X.stats$center, `-`, check.margin	= FALSE)
    smd <- sweep(smd, 2L, X.stats$scale, `/`, check.margin = FALSE)

    rowMeans(abs(smd))
  }

  lambda <- validate_lambda(lambda)
  W.hat <- predict(object, X.test, lambda = lambda, type = "response", .simplify = FALSE)

  loss0 <- loss1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    loss0 <- .imbalance(1 - W.test, 1 - W.hat$control)
  }
  if (!is.null(object[["_fit"]]$treated)) {
    loss1 <- .imbalance(W.test, W.hat$treated)
  }
  out <- list(control = loss0, treated = loss1)

  out[!vapply(out, is.null, logical(1))]
}
