## Pathwise Estimation of Covariate Balancing Propensity Scores

🚧 Work in progress 🚧

A package for pathwise estimation of regularized logistic propensity score models using covariate balancing loss functions rather than maximum likelihood. Regularization paths are fit via the [`adelie`](https://jamesyang007.github.io/adelie/) elastic-net solver with an interface inspired by [`glmnet`](https://glmnet.stanford.edu/), and objectives that directly target covariate balance for the ATE and ATT.

Some helpful links for getting started:

* [An introduction to `balnet`](https://erikcs.github.io/balnet/get-started.html)
* [Package docs](https://erikcs.github.io/balnet/reference.html)

### Installation

The development version can be installed via

```R
devtools::install_github("erikcs/balnet", subdir = "r-package/balnet")
```

Installing from source requires a C++17 compiler or later. To build with multithreading enabled, OpenMP needs to be available (on Mac, a simple option is to set the default C++ compiler to `gcc` installed via `brew`).
