# Simulated trade network (longitudinal)

A named list of three 10x10 matrices representing simulated bilateral
trade intensity between 10 countries across three time periods (2015,
2017, 2019). Useful for demonstrating longitudinal SRM analysis.

## Usage

``` r
data(trade_net)
```

## Format

A named list of three 10x10 numeric matrices. Countries: USA, CHN, DEU,
JPN, GBR, FRA, KOR, IND, BRA, CAN.

## Examples

``` r
data(trade_net)
fit = srm(trade_net)
summary(fit)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#> Actor                      0.5879    33.0%
#> Partner                    0.5259    29.5%
#> Unique                     0.6663    37.4%
#> Relationship (cov)         0.3575       --
#> Actor-Partner (cov)        0.2259       --
#> 
#> (Averaged across 3 time points)
```
