# balnet

Pathwise estimation of covariate balancing propensity scores using the [adelie](https://jamesyang007.github.io/adelie/) lasso solver.

> 🚧 Work in progress

## Installation

The development version can be installed via

```R
devtools::install_github("erikcs/balnet", subdir = "r-package/balnet")
```

Installing from source requires a C++17 compiler or later. To build with multithreading enabled, OpenMP needs to be available (on Mac, a simple option is compiling with `gcc` installed via `brew`).
