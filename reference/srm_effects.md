# Calculate SRM actor, partner, or unique effects

Decomposes a round-robin network into actor (sender), partner
(receiver), and unique (dyadic) effects using the Social Relations Model
framework. Requires a square (unipartite) adjacency matrix with at least
3 actors; matrices with fewer than 3 actors return `NA` values. For
bipartite networks, use
[`srm()`](https://netify-dev.github.io/srm/reference/srm.md) instead.

## Usage

``` r
srm_effects(mat, type = c("actor", "partner", "unique"), time = NULL)
```

## Arguments

- mat:

  A square numeric matrix (or list of square matrices for longitudinal
  data). Diagonal values are set to zero internally. NAs are replaced
  with zero. Must not be an `srm` object (effects are already stored in
  `fit$actor_effects`, etc.).

- type:

  One of `"actor"`, `"partner"`, or `"unique"`:

  - `"actor"`: sender effects, returned as an n x 1 matrix.

  - `"partner"`: receiver effects, returned as an n x 1 matrix.

  - `"unique"`: dyadic residuals, returned as an n x n matrix with zero
    diagonal.

- time:

  Optional character vector of time-point names to subset when `mat` is
  a named list.

## Value

For a single matrix, an effect matrix (with class `actor_effect`,
`partner_effect`, or `unique_effect`). For a list of matrices, a list of
such matrices. See
[`srm()`](https://netify-dev.github.io/srm/reference/srm.md) for the
unified interface.

## Details

When the input matrix is symmetric (undirected), actor and partner
effects are identical by construction. The decomposition is most
informative for directed (asymmetric) networks where sender and receiver
roles are distinct.

## Examples

``` r
# simple inline example
test_matrix = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow = 3, ncol = 3)
rownames(test_matrix) = colnames(test_matrix) = c("A", "B", "C")
actor_effect = srm_effects(test_matrix, type = "actor")

# cross-sectional with included dataset
data(classroom)
actor_eff = srm_effects(classroom, type = "actor")
head(actor_eff)
#>          actor_effect
#> Cady        1.0130833
#> Aaron       0.2873333
#> Gretchen   -0.2332500
#> Karen       1.0540833
#> Janis      -0.4112500
#> Damian      2.3567500

partner_eff = srm_effects(classroom, type = "partner")
unique_eff = srm_effects(classroom, type = "unique")

# longitudinal
data(trade_net)
actor_long = srm_effects(trade_net, type = "actor")
actor_2017 = srm_effects(trade_net, type = "actor", time = "2017")
```
