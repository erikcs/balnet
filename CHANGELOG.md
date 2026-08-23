# Changelog
All notable changes to `balnet` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.4] - 2026-08-??

### Fixed
- `cv.balnet`: Fix fold aggregation bug, fix folds argument handling, and add type.measure = "imbalance". [#45](https://github.com/erikcs/balnet/pull/45), [@ee36f54](https://github.com/erikcs/balnet/commit/ee36f54adac7a32fa32eec5bfe8c56782d06b247), [#48](https://github.com/erikcs/balnet/pull/48)
- Fix ridge lambda sequence discontinuity and increase ridge_scale to 0.05. [#54](https://github.com/erikcs/balnet/pull/54)
- Fix upstream adelie solver patch in 1-ULP/kkt screen. [#50](https://github.com/erikcs/balnet/pull/50)
- Fix mu0 prediction for linear predictor type. [@308c974](https://github.com/erikcs/balnet/commit/308c9747cf8c698a9fe05e05773b8c9946c79976)


## [0.0.3] - 2026-05-26

### Fixed
- Updated replication materials and improved print and summary methods.

## [0.0.2] - 2026-05-04

### Fixed
- Fix lambda sequence calculation when `alpha = 0`. [#38](https://github.com/erikcs/balnet/pull/38)

## [0.0.1] - 2026-04-03

Initial CRAN beta release.

## [0.0.0] - 2026-02-23

Beta version (source install only) released.
