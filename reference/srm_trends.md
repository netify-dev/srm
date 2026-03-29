# Track SRM effects over time

Extracts actor or partner effects across time points from a fitted `srm`
object and returns a tidy data frame suitable for plotting temporal
trends.

## Usage

``` r
srm_trends(object, type = c("actor", "partner"), actors = NULL)
```

## Arguments

- object:

  An object of class `srm` with multiple time points.

- type:

  One of `"actor"` or `"partner"`.

- actors:

  Optional character vector of actor names to include. If `NULL`, all
  actors are included.

## Value

A data frame with columns: `actor`, `time`, `effect`.

## Examples

``` r
# \donttest{
sim = sim_srm(n_actors = 6, n_time = 5, seed = 6886)
fit = srm(sim$Y)
trends = srm_trends(fit, type = "actor")
head(trends)
#>   actor time      effect
#> 1   n01   t1  1.09381997
#> 2   n02   t1  1.49088028
#> 3   n03   t1 -1.83100076
#> 4   n04   t1  0.12870414
#> 5   n05   t1 -0.81020547
#> 6   n06   t1 -0.07219816
# }
```
