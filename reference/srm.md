# Fit a Social Relations Model to network data

Main entry point for the `srm` package. Decomposes a network (or list of
networks) into actor, partner, and unique dyadic effects and computes
the full set of SRM variance and covariance components.

## Usage

``` r
srm(mat, time = NULL)
```

## Arguments

- mat:

  A square matrix, a list of square matrices (longitudinal), a
  rectangular matrix (bipartite), or a `netify` object. Diagonal values
  are set to zero and NAs are replaced with zero internally.

- time:

  Optional character vector of time-point names to subset to (only
  meaningful when `mat` is a named list).

## Value

An object of class `srm`, which is a list containing:

- actor_effects:

  Actor (sender) effects for each time point.

- partner_effects:

  Partner (receiver) effects for each time point.

- unique_effects:

  Unique dyadic effects for each time point.

- stats:

  A data frame of variance / covariance components for each time point.

- grand_mean:

  The overall network mean for each time point.

- matrices:

  The (possibly time-filtered) input matrices used.

- n_time:

  Number of time points.

- n_actors:

  Number of actors (per time point).

- bipartite:

  Logical; whether the model was fit as bipartite.

## Details

The SRM decomposes each observation \\X\_{ij}\\ as: \$\$X\_{ij} = \mu +
a_i + b_j + g\_{ij}\$\$ where \\\mu\\ is the grand mean, \\a_i\\ is the
actor (sender) effect, \\b_j\\ is the partner (receiver) effect, and
\\g\_{ij}\\ is the unique dyadic effect.

For unipartite (square) networks, five variance/covariance components
are estimated: actor variance, partner variance, unique variance,
relationship covariance (generalized reciprocity), and actor-partner
covariance. These require at least 4 actors (n \>= 4). With n = 3,
effects are estimated but variance components are returned as `NA`.

Because the SRM uses method-of-moments estimation, variance estimates
are not constrained to be non-negative. Negative estimates can occur
with small networks (roughly n \< 15) and indicate that the component
cannot be estimated reliably. When computing variance partitions,
negative estimates are clipped to zero before calculating percentages.

For bipartite (rectangular) networks, only actor variance, partner
variance, and unique variance are estimated (covariance components are
not defined when senders and receivers come from different populations).
Effects are simple deviations from the grand mean. Variance components
are bias-corrected using two-way ANOVA degrees of freedom. There is no
minimum size restriction beyond n \>= 2 in each mode.

**Symmetric networks:** When the input matrix is symmetric (undirected),
actor and partner effects are identical by construction. The variance
decomposition is constrained accordingly. The SRM is most informative
for directed (asymmetric) networks.

## See also

[`srm_effects()`](https://netify-dev.github.io/srm/reference/srm_effects.md)
for effect extraction,
[`srm_stats()`](https://netify-dev.github.io/srm/reference/srm_stats.md)
for individual statistics,
[`permute_srm()`](https://netify-dev.github.io/srm/reference/permute_srm.md)
for inference,
[`srm_trends()`](https://netify-dev.github.io/srm/reference/srm_trends.md)
and
[`srm_stability()`](https://netify-dev.github.io/srm/reference/srm_stability.md)
for longitudinal analysis,
[`sim_srm()`](https://netify-dev.github.io/srm/reference/sim_srm.md) for
simulation.

## Examples

``` r
mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
rownames(mat) = colnames(mat) = c("A", "B", "C")
fit = srm(mat)
fit
#> Social Relations Model
#> ---------------------------------------- 
#> Mode:       unipartite
#> Actors:     3
#> Grand mean: 1.5000
#> 
#> Variance Decomposition:
#>   (Not estimable -- network too small, need n >= 4)

data(classroom)
fit = srm(classroom)
fit
#> Social Relations Model
#> ---------------------------------------- 
#> Mode:       unipartite
#> Actors:     12
#> Grand mean: 3.5077
#> 
#> Variance Decomposition:
#>   Actor          1.4000  ( 42.3%)
#>   Partner        0.8539  ( 25.8%)
#>   Unique         1.0546  ( 31.9%)
#>   Relationship   0.6885  (cov)
#>   Actor-Partner  -0.2525  (cov)
summary(fit)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#> Actor                      1.4000    42.3%
#> Partner                    0.8539    25.8%
#> Unique                     1.0546    31.9%
#> Relationship (cov)         0.6885       --
#> Actor-Partner (cov)       -0.2525       --
```
