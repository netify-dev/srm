# Simulated classroom friendship ratings

A 12x12 matrix of simulated friendship ratings between students at North
Shore High, inspired by *Mean Girls*. Generated with actor, partner, and
unique effect components to serve as a simple cross-sectional example
with directed (asymmetric) ties.

## Usage

``` r
data(classroom)
```

## Format

A 12x12 numeric matrix with row and column names corresponding to
character names: Cady, Aaron, Gretchen, Karen, Janis, Damian, Kevin,
Glen, Shane, Trang, Regina, Norbury. Diagonal is zero.

## Examples

``` r
data(classroom)
fit = srm(classroom)
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
