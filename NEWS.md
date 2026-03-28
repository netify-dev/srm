# srm 1.0.0

## New Features

* **S3 class system**: New `srm()` function serves as the main entry point,
  returning an object of class `srm` with `print()`, `summary()`, and
  `plot()` methods.

* **Variance decomposition**: `summary.srm()` produces a full variance
  partition table showing actor, partner, and unique variance components
  with percentages.

* **Variance partition plot**: `plot(fit, type = "variance")` creates a
  bar chart of variance proportions.

* **Simulation**: `sim_srm()` generates synthetic network data from known
  SRM parameters for validation, teaching, and power analysis. Supports
  cross-sectional, longitudinal, and bipartite designs.

* **Permutation inference**: `permute_srm()` provides permutation-based
  tests for whether variance components significantly differ from zero,
  with `print()` and `plot()` methods for the results.

* **Bipartite (two-mode) networks**: All functions now support rectangular
  (bipartite) adjacency matrices where senders and receivers come from
  different populations.

* **Longitudinal analysis**:
    - `srm_trends()` extracts actor/partner effects over time as a
      tidy data frame.
    - `srm_trend_plot()` visualizes temporal trajectories.
    - `srm_stability()` computes rank-order correlations between
      consecutive time points.

* **Improved netify integration**: `srm()` directly accepts `netify`
  objects (both cross-sectional and longitudinal) without manual
  conversion.

* **Additional datasets**:
    - `classroom`: Simulated 12-student friendship ratings.
    - `trade_net`: Simulated 10-country trade network over 3 periods.

## Documentation

* Added methodology vignette with full mathematical framework.
* Added pipeline vignette showing netify-srm integration.
* Added pkgdown website at https://netify-dev.github.io/srm/.
* Added GitHub Actions for CI/CD (pkgdown deployment, release binaries).

## Infrastructure

* Package version bumped to 0.1.0.
* Added URL and BugReports fields to DESCRIPTION.
* Expanded test suite from 3 to 12 test files.
* Added `stats` to Imports.
