#' Compute weighted column means and standard deviations.
#'
#' For X column `x`, and weight column `w`, this function computes
#' \eqn{\mu = \frac{\sum_{i}^{n} w_i x_i}{\sum_{i}^{n} w_i}} and
#' \eqn{\sigma^2 = \frac{\sum_{i}^{n} w_i (x_i - \mu)^2}{\sum_{i}^{n} w_i}.}
#'
#' @param X A `n * p` numeric R matrix.
#' @param weights A `n * L` numeric weight matrix, or an `n`-vector.
#' @param compute_sd Whether to return the standard deviation.
#'
#' @return `L * p` matrices of column stats.
#'
#' @keywords internal
col_stats <- function(
  X,
  weights = NULL,
  compute_sd = FALSE,
  n_threads = 1L
)
{
  if (is.null(weights)) {
    weights <- matrix(1, nrow(X), 1L)
  } else {
    weights <- as.matrix(weights)
  }
  stopifnot(nrow(X) == nrow(weights))

  rcpp_col_stats(X, weights, compute_sd, n_threads)
}

#' Quickly center and scale a standard dense R matrix.
#'
#' @param X A numeric R matrix.
#' @param weights Sample weights.
#' @param standardize Whether to standardize.
#' @param inplace Whether to overwrite X.
#' @param n_threads Number of threads used.
#'
#' @return A list containing the standardized information.
#'
#' @keywords internal
standardize <- function(
  X,
  weights = NULL,
  standardize = TRUE,
  inplace = FALSE,
  n_threads = 1L
)
{
  if (standardize) {
    if (is.null(weights)) {
      weights <- matrix(1, nrow(X), 1L)
    } else {
      weights <- as.matrix(weights)
    }
    col.stats <- rcpp_col_stats(X, weights, compute_sd = TRUE, n_threads = n_threads)
    center <- drop(col.stats$center)
    scale <- drop(col.stats$scale)
    scale[scale <= 0] <- 1

    if (inplace) {
      rcpp_standardize_inplace(X, center, scale, n_threads)
    } else {
      X <- rcpp_standardize(X, center, scale, n_threads)
    }
  } else {
    center <- rep(0.0, ncol(X))
    scale <- rep(1.0, ncol(X))
  }

  list(X = X, center = center, scale = scale)
}

sp_tcrossprod <- function(
  X,
  beta,
  n_threads = 1L
) {
  stopifnot(ncol(X) == ncol(beta))
  stopifnot(is.matrix(X))
  stopifnot(is(beta, "dgRMatrix"))

  rcpp_sp_tcrossprod(X, beta, n_threads)
}

get_lambda_min_ratio <- function(lambda.min.ratio, max.imbalance, X.stan, W, sample.weights, target, alpha) {
  if (is.null(max.imbalance)) {
    out <- c(lambda.min.ratio, lambda.min.ratio)
  } else {
    if (alpha < 1) {
      stop("Setting max.imbalance is only possible with lasso (alpha = 1).")
    }
    if (max.imbalance <= 0) {
      stop("max.imbalance should be > 0.")
    }
    lambda.min.ratio0 <- lambda.min.ratio1 <- lambda.min.ratio
    if (target %in% c("ATE", "ATT", "control")) {
      stats0 <- col_stats(X.stan, weights = (1 - W) * sample.weights)
      lambda0.max <- max(abs(stats0$center)) # Note, this assumes X.stan is standardized.
      if (max.imbalance < lambda0.max) {
        lambda.min.ratio0 <- max.imbalance / lambda0.max
      }
    }
    if (target %in% c("ATE", "treated")) {
      stats1 <- col_stats(X.stan, weights = W * sample.weights)
      lambda1.max <- max(abs(stats1$center))
      if (max.imbalance < lambda1.max) {
        lambda.min.ratio1 <- max.imbalance / lambda1.max
      }
    }
    out <- c(lambda.min.ratio0, lambda.min.ratio1)
  }

  out
}

validate_lambda <- function(lambda) {
  if (is.character(lambda)) {
    stop("Unsupported lambda argument.")
  }

  if (is.null(lambda)) {
    lambda.in <- list(NULL, NULL)
  } else if (is.list(lambda)) {
    if (length(lambda) < 2) {
     lambda.in <- list(lambda[[1]], lambda[[1]])
    } else {
      lambda.in <- lambda
    }
  } else if (is.matrix(lambda)) {
    if (ncol(lambda) < 2) {
      lambda.in <- list(lambda[, 1], lambda[, 1])
    } else {
      lambda.in <- list(lambda[, 1], lambda[, 2])
    }
  } else if (is.vector(lambda)) {
    lambda.in <- list(lambda, lambda)
  } else {
    stop("Unsupported lambda argument.")
  }

  lambda.in
}

validate_groups <- function(groups, p, colnames) {
  if (is.null(groups)) return(NULL)

  if (!is.list(groups) || length(groups) == 0) {
    stop("groups must be a list of integer vectors.")
  }

  all.indices <- unlist(groups)
  if (!is.numeric(all.indices) || any(all.indices < 1) || any(all.indices > p)) {
    stop("group indices must be between 1 and ncol(X).")
  }
  if (anyDuplicated(all.indices)) {
    stop("Indices cannot belong to multiple groups.")
  }

  for (i in seq_along(groups)) {
    g <- groups[[i]]
    if (length(g) > 1) {
      if (is.unsorted(g)) {
        stop("group indices need to be in increasing order.")
      }
      if (any(diff(g) != 1)) {
        stop("groups should be contiguous column indices.")
      }
    }
  }

  if (any(names(groups) %in% colnames)) {
    stop("group names cannot equal individual column names of X.")
  }
}

process_groups <- function(groups, p) {
  if (is.null(groups)) {
    groups <- 0:(p - 1)
    group_sizes <- rep_len(1L, length(groups))
  } else {
    grouped.starts <- sapply(groups, min)
    solo.cols <- setdiff(1:p, unlist(groups))

    starts.1indexed <- sort(c(grouped.starts, solo.cols))
    temp.bounds <- c(starts.1indexed, p + 1)

    groups <- as.integer(starts.1indexed - 1)
    group_sizes <- as.integer(diff(temp.bounds))
  }

  list(
    groups = groups,
    group_sizes = group_sizes,
    G = length(groups)
  )
}

# Construct a p * p.new sparse aggregation matrix mapping original columns
# to group-level means (weights 1/|group|), with identity mapping
# for columns not assigned to any group.
col_group_mat <- function(groups, p, colnames) {
  validate_groups(groups, p, colnames)

  all.idx <- unlist(groups)
  solo.indices <- setdiff(1:p, all.idx)

  group.defs <- lapply(seq_along(groups), function(k) {
    idx <- groups[[k]]
    nm <- names(groups)[k]
    if (is.null(nm) || nm == "") {
     nm <- paste0("Grp_", idx[1], "-", idx[length(idx)])
    }
    list(idx = idx, name = nm, sort.key = min(idx))
  })

  solo.defs <- lapply(solo.indices, function(i) {
    list(idx = i, name = colnames[i], sort.key = i)
  })

  # Merge and sort (interleave groups into their natural position)
  # sort by where the feature first appears in the original matrix
  all.defs <- c(group.defs, solo.defs)
  all.defs <- all.defs[order(sapply(all.defs, function(x) x$sort.key))]

  p.new <- length(all.defs)
  i.list <- vector("list", p.new)
  j.list <- vector("list", p.new)
  x.list <- vector("list", p.new)
  new.names <- character(p.new)
  for (col.new in seq_len(p.new)) {
    def <- all.defs[[col.new]]
    rows <- def$idx
    n.rows <- length(rows)

    i.list[[col.new]] <- rows
    j.list[[col.new]] <- rep(col.new, n.rows)
    x.list[[col.new]] <- rep(1/n.rows, n.rows)
    new.names[col.new] <- def$name
  }

  Matrix::sparseMatrix(
    i = unlist(i.list),
    j = unlist(j.list),
    x = unlist(x.list),
    dims = c(p, p.new),
    dimnames = list(NULL, new.names)
  )
}

collapse_X <- function(X, groups, colnames) {
  M <- col_group_mat(groups, ncol(X), colnames)

  as.matrix(X %*% M)
}

# From glmnet: linear interpolation of lambda values
lambda.interp = function(lambda, s) {
  ### lambda is the index sequence that is produced by the model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indices, and a vector of fractions.
  ### the new values are interpolated bewteen the two using the fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)

  if (length(lambda) == 1) { # degenerate case of only one lambda
    nums = length(s)
    left = rep(1, nums)
    right = left
    sfrac = rep(1, nums)
  }
  else {
    ## s[s > max(lambda)] = max(lambda)
    ## s[s < min(lambda)] = min(lambda)
    k = length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    coord <- stats::approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac = (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] = 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] = 1

  }
  list(left = left, right = right, frac = sfrac)
}
