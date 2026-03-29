# Simulate network data from a Social Relations Model

Generates synthetic network matrices from known actor, partner, and
unique effects, useful for teaching, validation, and power analysis.

## Usage

``` r
sim_srm(
  n_actors = 10,
  n_time = 1,
  actor_var = 1,
  partner_var = 1,
  unique_var = 1,
  actor_partner_cov = 0,
  relationship_cov = 0,
  grand_mean = 0,
  bipartite = FALSE,
  seed = NULL
)
```

## Arguments

- n_actors:

  Number of actors (nodes) in the network. For bipartite networks, a
  length-2 vector `c(n_row, n_col)` specifying the number of senders and
  receivers.

- n_time:

  Number of time periods (default 1 for cross-sectional).

- actor_var:

  Variance of actor (sender) effects.

- partner_var:

  Variance of partner (receiver) effects.

- unique_var:

  Variance of unique dyadic effects.

- actor_partner_cov:

  Covariance between actor and partner effects. Must satisfy
  `actor_partner_cov^2 <= actor_var * partner_var`.

- relationship_cov:

  Covariance of unique effects within dyads (generalized reciprocity).
  Must satisfy `abs(relationship_cov) <= unique_var`.

- grand_mean:

  Overall network mean.

- bipartite:

  If `TRUE`, generate a bipartite (two-mode) network. When bipartite,
  `n_actors` is a length-2 vector `c(n_row, n_col)`. The
  `actor_partner_cov` and `relationship_cov` parameters are ignored for
  bipartite networks (these quantities are not defined when senders and
  receivers are from different populations).

- seed:

  Optional random seed for reproducibility.

## Value

A list with components:

- Y:

  A matrix or named list of matrices (if `n_time > 1`).

- truth:

  A list containing the true `actor_effects`, `partner_effects`,
  `unique_effects`, and `grand_mean`.

- params:

  The parameter values used for simulation.

## Examples

``` r
# cross-sectional
sim = sim_srm(n_actors = 10, actor_var = 1, partner_var = 0.5,
               unique_var = 2, seed = 6886)
fit = srm(sim$Y)
summary(fit)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#> Actor                      0.4843    15.2%
#> Partner                    0.6512    20.4%
#> Unique                     2.0533    64.4%
#> Relationship (cov)        -0.1958       --
#> Actor-Partner (cov)        0.5433       --

# \donttest{
# longitudinal
sim_long = sim_srm(n_actors = 8, n_time = 5, seed = 6886)
fit_long = srm(sim_long$Y)
fit_long
#> Social Relations Model
#> ---------------------------------------- 
#> Mode:       unipartite
#> Time points: 5
#> Actors:      8
#> Grand mean:  -0.4985 (avg)
#> 
#> Variance Decomposition:
#>   Actor          0.8418  ( 29.1%)
#>   Partner        1.0098  ( 34.9%)
#>   Unique         1.0388  ( 35.9%)
#>   Relationship  -0.0524  (cov)
#>   Actor-Partner   0.0750  (cov)

# bipartite
sim_bip = sim_srm(n_actors = c(6, 8), bipartite = TRUE, seed = 6886)
fit_bip = srm(sim_bip$Y)
summary(fit_bip)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#> Actor                      1.2467    39.3%
#> Partner                    0.8214    25.9%
#> Unique                     1.1075    34.9%
# }
```
