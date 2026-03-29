# Plot SRM effects over time

Creates a line plot of actor or partner effects across time points.

## Usage

``` r
srm_trend_plot(object, type = c("actor", "partner"), actors = NULL, n = 8)
```

## Arguments

- object:

  An object of class `srm` with multiple time points.

- type:

  One of `"actor"` or `"partner"`.

- actors:

  Optional character vector of actor names to highlight. If `NULL`, all
  actors are shown.

- n:

  If `actors` is `NULL`, the top `n` actors by mean absolute effect are
  shown (default 8).

## Value

A `ggplot` object.

## Examples

``` r
# \donttest{
sim = sim_srm(n_actors = 8, n_time = 5, seed = 6886)
fit = srm(sim$Y)
srm_trend_plot(fit, type = "actor", n = 4)

# }
```
