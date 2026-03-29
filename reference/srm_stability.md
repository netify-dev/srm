# Compute stability of effects across time

Calculates correlations of actor (or partner) effects between
consecutive time points, measuring how stable individual positions are
over time.

## Usage

``` r
srm_stability(object, type = c("actor", "partner"))
```

## Arguments

- object:

  An object of class `srm` with multiple time points.

- type:

  One of `"actor"` or `"partner"`.

## Value

A data frame with columns: `time1`, `time2`, `correlation`, `n`. The
`correlation` column is `NA` when fewer than 3 actors are shared between
consecutive time points (Pearson correlation requires at least 3
observations).

## Examples

``` r
# \donttest{
sim = sim_srm(n_actors = 10, n_time = 4, seed = 6886)
fit = srm(sim$Y)
srm_stability(fit, type = "actor")
#>   time1 time2 correlation  n
#> 1    t1    t2  0.06041908 10
#> 2    t2    t3  0.07972254 10
#> 3    t3    t4  0.43008139 10
# }
```
