## Pathwise Estimation of Covariate Balancing Propensity Scores

A package for pathwise estimation of regularized logistic propensity score models using covariate balancing loss functions rather than maximum likelihood. Regularization paths are fit via the [adelie](https://jamesyang007.github.io/adelie/) elastic-net solver with a *glmnet*-like interface and objectives that directly target covariate balance for the ATE and ATT.

Some useful links:

* [An introduction to balnet](https://erikcs.github.io/balnet/get-started.html)
* [Package docs](https://erikcs.github.io/balnet/reference.html)

### Installation

The development version can be installed via

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
# Note: lambda = 0 selects the final lambda in the sequence. Scalar values
# are applied to both arms.
plot(fit, lambda = 0)

# Predict propensity scores at end of lambda path.
W.hat <- predict(fit, X, lambda = 0)

# Get weights at end of lambda path.
ipw.weights <- ipw(fit, lambda = 0)

# Estimate ATE using IPW weights.
mean(Y * (ipw.weights$treated - ipw.weights$control))
```

### References

Erik Sverdrup and Trevor Hastie.
<b>balnet: Pathwise Estimation of Covariate Balancing Propensity Scores.</b> 2026.
[<a href="https://arxiv.org/abs/2602.18577">arxiv</a>]
