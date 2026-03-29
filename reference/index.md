# Package index

## Model Fitting

Fit Social Relations Models to network data

- [`srm()`](https://netify-dev.github.io/srm/reference/srm.md) : Fit a
  Social Relations Model to network data
- [`srm_effects()`](https://netify-dev.github.io/srm/reference/srm_effects.md)
  : Calculate SRM actor, partner, or unique effects
- [`srm_stats()`](https://netify-dev.github.io/srm/reference/srm_stats.md)
  : Calculate SRM summary statistics

## S3 Methods

Standard R methods for srm objects

- [`print(`*`<srm>`*`)`](https://netify-dev.github.io/srm/reference/print.srm.md)
  : Print an srm object
- [`summary(`*`<srm>`*`)`](https://netify-dev.github.io/srm/reference/summary.srm.md)
  : Summarize an srm object
- [`print(`*`<summary.srm>`*`)`](https://netify-dev.github.io/srm/reference/print.summary.srm.md)
  : Print a summary.srm object
- [`plot(`*`<srm>`*`)`](https://netify-dev.github.io/srm/reference/plot.srm.md)
  : Plot an srm object

## Inference

Permutation-based significance testing

- [`permute_srm()`](https://netify-dev.github.io/srm/reference/permute_srm.md)
  : Permutation test for SRM variance components
- [`print(`*`<srm_permtest>`*`)`](https://netify-dev.github.io/srm/reference/print.srm_permtest.md)
  : Print a permutation test result
- [`plot(`*`<srm_permtest>`*`)`](https://netify-dev.github.io/srm/reference/plot.srm_permtest.md)
  : Plot permutation test null distributions

## Longitudinal Analysis

Track effects over time

- [`srm_trends()`](https://netify-dev.github.io/srm/reference/srm_trends.md)
  : Track SRM effects over time
- [`srm_trend_plot()`](https://netify-dev.github.io/srm/reference/srm_trend_plot.md)
  : Plot SRM effects over time
- [`srm_stability()`](https://netify-dev.github.io/srm/reference/srm_stability.md)
  : Compute stability of effects across time

## Visualization

Additional plotting functions

- [`srm_plot()`](https://netify-dev.github.io/srm/reference/srm_plot.md)
  : Create plots for SRM actor, dyadic, partner effects

## Simulation

Generate synthetic network data

- [`sim_srm()`](https://netify-dev.github.io/srm/reference/sim_srm.md) :
  Simulate network data from a Social Relations Model

## Datasets

Example datasets

- [`atop_EA`](https://netify-dev.github.io/srm/reference/atop_EA.md) :
  East Asian consultation agreements from ATOP
- [`classroom`](https://netify-dev.github.io/srm/reference/classroom.md)
  : Simulated classroom friendship ratings
- [`small_council`](https://netify-dev.github.io/srm/reference/small_council.md)
  : Simulated Small Council pursuit data (bipartite)
- [`trade_net`](https://netify-dev.github.io/srm/reference/trade_net.md)
  : Simulated trade network (longitudinal)
