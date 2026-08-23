#' Tuning via subsampling for balnet.
#'
#' Tunes a balnet object using sub sampled imbalances calculated using full-sample weights
#' a la Wang & Zubizarreta (2020).
#'
#' @param X A numeric matrix or data frame with pre-treatment covariates.
#' @param W Treatment vector (0: control, 1: treated).
#' @param B The number of bootstrap replicates.
#' @param ... Arguments for \code{\link{balnet}}.
#'
#' @return A fit cv.boot.balnet object (same functionality as cv.balnet)
#'
#' @references Wang, Yixin, and Jose R. Zubizarreta.
#'  "Minimal dispersion approximately balancing weights: asymptotic properties and practical considerations."
#'  Biometrika 107(1), 2020.
#'
#' @examples
#' \donttest{
#' n <- 100
#' p <- 25
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 1 / (1 + exp(1 - X[, 1])))
#'
#' # Fit an ATE model.
#' cv.fit <- cv.boot.balnet(X, W)
#'
#' # Print CV summary.
#' print(cv.fit)
#'
#' # Plot at tuned lambda.
#' plot(cv.fit)
#'
#' # Predict at tuned lambda.
#' W.hat <- predict(cv.fit, X)
#'
#' }
#'
#' @export
cv.boot.balnet <- function(
  X,
  W,
  B = 500,
  ...
)
{
  dot.args <- list(...)
  get_loss <- get_imbalance
  X.stats <- col_stats(X, dot.args[["sample.weights"]], compute_sd = TRUE)
  X.stats$scale[X.stats$scale <= 0] <- 1
  test.list <- replicate(B, sample.int(length(W), length(W) %/% 2), simplify = FALSE)


  if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]]) cat("Fitting full model\n")
  fit.full <- balnet(X, W, ...)
  lambda.full <- fit.full[["_lambda"]]
  sample.weights <- fit.full[["sample.weights"]]

  W.hat.full <- predict(fit.full, X, lambda = lambda.full, type = "response", .simplify = FALSE)

  cv.list <- list()
  for (k in 1:length(test.list)) {
    if (!is.null(dot.args[["verbose"]]) && dot.args[["verbose"]] && refit) cat(sprintf("\nFold: %d/%d\n", k, nfolds))
    test <- test.list[[k]]
    X.test <- X[test, , drop = FALSE]
    W.test <- W[test]
    sample.weights.test <- sample.weights[test]
    W.hat.test <- lapply(W.hat.full, function(m) m[test, , drop = FALSE])
    loss <- do.call(
        get_loss,
        list(fit.full, X.test, W.test, sample.weights.test, lambda.full, X.stats = X.stats, W.hat = W.hat.test)
    )

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
    "type.measure" = "bootstrapped.imbalance"
  )

  fit.full[["lambda.min"]] <- if (length(lambda.min.out) > 1) lambda.min.out else lambda.min.out[[1]]
  fit.full[["_cv.info"]] <- cv.info
  fit.full[["call"]] <- match.call()
  class(fit.full) <- c("cv.boot.balnet", "cv.balnet", class(fit.full))

  fit.full
}
