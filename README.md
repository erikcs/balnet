## balnet

[![CRANstatus](https://www.r-pkg.org/badges/version/balnet)](https://cran.r-project.org/package=balnet)

A package for pathwise estimation of regularized logistic propensity score models using covariate balancing loss functions rather than maximum likelihood. Regularization paths are fit via the [adelie](https://jamesyang007.github.io/adelie/) elastic-net solver with a *glmnet*-like interface, yielding balancing weights that target covariate balance for the ATE and ATT.

The choice of penalty determines which norm of the covariate imbalance vector is bounded by $\lambda$, and so the regularization path traces a sequence of decreasing imbalance tolerances. Lasso bounds the maximum covariate imbalance ($\ell_\infty$) and ridge bounds the $\ell_2$ norm.

Some useful links:

* [An introduction to balnet](https://erikcs.github.io/balnet/balnet.html)
* [Package manual](https://erikcs.github.io/balnet/reference.html)

### Installation

The latest release of the package can be installed through CRAN:

```R
install.packages("balnet")
```

The development version can be installed via:

```R
devtools::install_github("erikcs/balnet", subdir = "r-package/balnet")
```

Installing from source requires a C++17 compiler or later. To build with multithreading enabled, OpenMP needs to be available (on Mac, a simple option is to set the default C++ compiler to `gcc` installed via `brew`).

### Usage Example

```R
# Simulate data with confounding.
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 1 / (1.5 + exp(X[, 2] + X[, 3])))
Y <- W + 2 * log(1 + exp(X[, 1] + X[, 2] + X[, 3])) + rnorm(n)

# Fit model targeting the ATE = E[Y(1)] - E[Y(0)].
# Two logistic models are fit: one for treated, one for control.
fit <- balnet(X, W, target = "ATE")

# Print path summary.
print(fit)

# Visualize the path.
plot(fit)

# Plot the covariate imbalance at given lambda.
# lambda = 0 selects the final lambda in the sequences.
plot(fit, lambda = 0)

# Predict propensity scores at end of lambda path.
W.hat <- predict(fit, X, lambda = 0)

# Get balancing weights at end of lambda path.
ipw.weights <- balweights(fit, lambda = 0)

# Estimate ATE using balancing weights.
mean(Y * (ipw.weights$treated - ipw.weights$control))
```

### References

Jerome H. Friedman, Trevor Hastie, and Rob Tibshirani.
<b>Regularization Paths for Generalized Linear Models via Coordinate Descent.</b>
<i>Journal of Statistical Software</i>, 33, 2010.
[<a href="https://www.jstatsoft.org/article/view/v033i01">paper</a>]

Erik Sverdrup and Trevor Hastie.
<b>balnet: Pathwise Estimation of Covariate Balancing Propensity Scores.</b> 2026.
[<a href="https://arxiv.org/abs/2602.18577">arxiv</a>]

James Yang and Trevor Hastie.
<b>A Fast and Scalable Pathwise-Solver for Group Lasso and Elastic Net Penalized Regression via Block-Coordinate Descent.</b> 2024.
[<a href="https://arxiv.org/abs/2405.08631">arxiv</a>]
