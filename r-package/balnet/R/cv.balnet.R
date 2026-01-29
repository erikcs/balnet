get_balance_loss <- function(object, X, W, sample.weights, lambda) {
  .balance_loss <- function(W) {
    colSums(sample.weights * (W * exp(-eta) + (1 - W) * eta)) / sum(sample.weights)
  }

  lambda <- validate_lambda(lambda)
  loss0 <- loss1 <- NULL
  if (!is.null(object[["_fit"]]$control)) {
    eta <- predict(object[["_fit"]]$control, X, lambda = lambda[[1]], type = "link")
    loss0 <- .balance_loss(1 - W)
  }
  if (!is.null(object[["_fit"]]$treated)) {
    eta <- predict(object[["_fit"]]$treated, X, lambda = lambda[[2]], type = "link")
    loss1 <- .balance_loss(W)
  }
  out <- list(control = loss0, treated = loss1)

  out[!vapply(out, is.null, logical(1))]
}

#' Cross-validation for balnet.
#'
#' @param X A numeric matrix or data frame with pre-treatment covariates.
#' @param W Treatment vector (0: control, 1: treated).
#' @param type.measure The loss to minimize for cross-validation. Default is balance loss.
#' @param nfolds The number of folds used for cross-validation, default is 10.
#' @param foldid An optional `n`-vector specifying which fold 1 to `nfold` a sample belongs to.
#' If NULL, this defaults to `sample(rep(seq(nfolds), length.out = nrow(X)))`.
#' @param ... Arguments for \code{\link{balnet}}.
#'
#' @return A fit cv.balnet object.
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
#' pp <- predict(cv.fit, X)
#'
#' # Extract coefficients at cross-validated lambda.
#' coefs <- coef(cv.fit)
#' }
#'
#' @export
cv.balnet <- function(
  X,
  W,
  type.measure = c("balance.loss"),
  nfolds = 10,
  foldid = NULL,
  ...
)
{
  type.measure <- match.arg(type.measure)
  if (type.measure == "balance.loss") {
    get_loss <- `get_balance_loss`
  }
  nfolds <- max(nfolds, 3)
  if (is.null(foldid)) {
    foldid <- sample(rep(seq(nfolds), length.out = nrow(X)))
  }
  dot.args <- list(...)

  if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) message("Fitting full model")
  fit.full <- balnet(X, W, ...)
  lambda.full <- fit.full[["lambda"]]
  sample.weights <- fit.full[["sample.weights"]]

  cv.list <- list()
  for (k in 1:nfolds) {
    if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) message(sprintf("\nFold: %d/%d", k, nfolds))
    test <- foldid == k
    train <- !test
    X.train <- X[train, , drop = FALSE]
    W.train <- W[train]
    fit.train <- balnet(X.train, W.train, standardize = "inplace", ...)

    X.test <- X[test, , drop = FALSE]
    W.test <- W[test]
    sample.weights.test <- sample.weights[test]
    loss <- do.call(get_loss, list(fit.train, X.test, W.test, sample.weights.test, lambda.full)) # TODO-balnet gradient norm loss
    cv.list[[k]] <- loss
  }
  cv.mean0 <- cv.mean1 <- NULL
  idx.min0 <- idx.min1 <- NA
  lambda.min0 <- lambda.min1 <- NA
  if (!is.null(cv.list[[1]][["control"]])) {
    cv.mean0 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, "control")), nfolds, length(lambda.full$control)))
    idx.min0 <- which.min(cv.mean0)
    lambda.min0 <- lambda.full[["control"]][idx.min0]
  }
  if (!is.null(cv.list[[1]][["treated"]])) {
    cv.mean1 <- colMeans(matrix(unlist(lapply(cv.list, `[[`, "treated")), nfolds, length(lambda.full$treated)))
    idx.min1 <- which.min(cv.mean1)
    lambda.min1 <- lambda.full[["treated"]][idx.min1]
  }
  cv.info <- list(
    "cv.mean" = list(control = cv.mean0, treated = cv.mean1),
    "idx.min" = list(control = idx.min0, treated = idx.min1),
    "lambda.min" = list(control = lambda.min0, treated = lambda.min1),
    "type.measure" = type.measure
  )

  fit.full[["cv.info"]] <- cv.info #TODO-balnet: store "keep"-like entry like glmnet?
  fit.full[["call"]] <- match.call()
  class(fit.full) <- c("cv.balnet", class(fit.full))

  fit.full
}

#' @rdname lambda
#' @method lambda balnet
#' @export
lambda.cv.balnet <- function(
  object,
  lambda = "lambda.min",
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- object[["cv.info"]]$lambda.min
  } else if (is.null(lambda)) {
    return(coef.balnet(object, lambda = lambda))
  } else {
    stop("Invalid lambda.")
  }
  out <- object[["lambda"]]
  out.nn <- out[!vapply(out, is.null, logical(1))]

  if (length(out.nn) > 1) {
    return(out.nn)
  } else {
    return(out.nn[[1]])
  }
}

#' Extract coefficients from a cv.balnet object.
#'
#' @param object A `cv.balnet` object.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
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
#' cv.fit <- cv.balnet(X, W)
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
    lambda <- object[["cv.info"]]$lambda.min
  }

  coef.balnet(object, lambda = lambda)
}

#' Predict using a cv.balnet object.
#'
#' @param object A `cv.balnet` object.
#' @param newx A numeric matrix.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
#' @param type The type of predictions. Default is "response" (propensity scores).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Estimated predictions. For dual-arm fits (control and treatment),
#'   returns a list containing predictions for each arm.
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
#' # Predict at cross-validated lambda.
#' pp <- predict(cv.fit, X)
#' }
#'
#' @method predict cv.balnet
#' @export
predict.cv.balnet <- function(
  object,
  newx,
  lambda = "lambda.min",
  type = c("response"),
  ...
)
{
  if (identical(lambda, "lambda.min")) {
    lambda <- object[["cv.info"]]$lambda.min
  }

  predict.balnet(object, newx, lambda = lambda, type = type)
}

#' Print a cv.balnet object.
#'
#' @param x A `cv.balnet` object.
#' @param digits Number of digits to print.
#' @param ... Additional print arguments.
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

  utils::capture.output(out <- print.balnet(x, digits = digits, drop = FALSE, ...))
  df0 <- df1 <- data.frame()
  if (!is.null(x[["_fit"]]$control)) {
    idx.min0 <- x[["cv.info"]]$idx.min[[1]]
    df0 <- cbind(Arm = "Control", out$control[idx.min0, ], Index = idx.min0)
  }
  if (!is.null(x[["_fit"]]$treated)) {
    idx.min1 <- x[["cv.info"]]$idx.min[[2]]
    df1 <- cbind(Arm = "Treated", out$treated[idx.min1, ], Index = idx.min1)
  }

  type.measure.nice <- gsub(".", " ", x[["cv.info"]]$type.measure, fixed = TRUE)
  cat("Lambda min (", type.measure.nice, "):\n", sep = "")
  print(rbind(df0, df1), digits = digits, row.names = FALSE, right = FALSE)
}

#' Plot diagnostics for a `cv.balnet` object.
#'
#' @param x A `cv.balnet` object.
#' @param lambda The lambda to use. Defaults to the cross-validated lambda.
#' @param ... Additional arguments.
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
    lambda <- x[["cv.info"]]$lambda.min
  }

  plot.balnet(x, lambda = lambda, ...)
}
