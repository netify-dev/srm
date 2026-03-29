# Calculate SRM summary statistics

Computes descriptive statistics and variance/covariance components from
a round-robin network using the Social Relations Model framework.
Requires a square (unipartite) adjacency matrix; for bipartite networks,
use [`srm()`](https://netify-dev.github.io/srm/reference/srm.md)
instead.

## Usage

``` r
srm_stats(
  mat,
  type = c("rowmeans", "colmeans", "totalmeans", "actor_var", "unique_var",
    "partner_var", "relationship_cov", "actor_partner_cov"),
  time = NULL
)
```

## Arguments

- mat:

  A square numeric matrix (or list of square matrices for longitudinal
  data). Diagonal values are set to zero internally. NAs are replaced
  with zero. Must not be an `srm` object (statistics are already in
  `summary(fit)$variance_table`).

- type:

  One of:

  - `"rowmeans"`: row means (average out-ties per actor).

  - `"colmeans"`: column means (average in-ties per actor).

  - `"totalmeans"`: grand mean of off-diagonal entries.

  - `"actor_var"`: actor (sender) variance component.

  - `"partner_var"`: partner (receiver) variance component.

  - `"unique_var"`: unique dyadic variance component.

  - `"relationship_cov"`: generalized reciprocity covariance.

  - `"actor_partner_cov"`: covariance between sending and receiving.

  Variance and covariance components require n \>= 4 actors. For n \<=
  3, `NA` is returned.

- time:

  Optional character vector of time-point names to subset when `mat` is
  a named list.

## Value

A named numeric vector (for descriptive statistics) or scalar (for
variance components). For a list input, a list of such values. See
[`srm()`](https://netify-dev.github.io/srm/reference/srm.md) for the
unified interface.

## Examples

``` r
# simple inline example
test_matrix = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow = 3, ncol = 3)
rownames(test_matrix) = colnames(test_matrix) = c("A", "B", "C")
actor_totalmeans = srm_stats(test_matrix, type = "totalmeans")

data(classroom)
srm_stats(classroom, type = "rowmeans")
#>      Cady     Aaron  Gretchen     Karen     Janis    Damian     Kevin      Glen 
#> 4.4172727 3.6690909 3.4300000 4.5300000 3.1218182 5.8963636 3.7109091 2.9563636 
#>     Shane     Trang    Regina   Norbury 
#> 2.6563636 3.0945455 0.7109091 3.8990909 
srm_stats(classroom, type = "actor_var")
#> [1] 1.400046
srm_stats(classroom, type = "actor_partner_cov")
#> [1] -0.2524702

# longitudinal
data(trade_net)
srm_stats(trade_net, type = "colmeans", time = c("2015", "2019"))
#> $`2015`
#>       USA       CHN       DEU       JPN       GBR       FRA       KOR       IND 
#> 2.2188889 0.8955556 2.3422222 2.6511111 3.0111111 1.7422222 1.9833333 2.4411111 
#>       BRA       CAN 
#> 1.6144444 3.6844444 
#> 
#> $`2019`
#>      USA      CHN      DEU      JPN      GBR      FRA      KOR      IND 
#> 2.585556 2.523333 2.006667 2.728889 1.484444 3.736667 1.736667 2.053333 
#>      BRA      CAN 
#> 2.988889 3.095556 
#> 
```
