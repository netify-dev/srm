# Summarize an srm object

Produces a detailed variance decomposition table and, for longitudinal
models, per-period breakdowns.

## Usage

``` r
# S3 method for class 'srm'
summary(object, ...)
```

## Arguments

- object:

  An object of class `srm`.

- ...:

  Additional arguments (ignored).

## Value

An object of class `summary.srm`, which is a list containing:

- stats:

  The full stats data frame.

- variance_table:

  A summary table with variance, percentage, and component labels.

- n_time:

  Number of time points.

- bipartite:

  Logical; bipartite indicator.

## Examples

``` r
mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
rownames(mat) = colnames(mat) = c("A", "B", "C")
fit = srm(mat)
summary(fit)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#>   (Not estimable -- network too small, need n >= 4)
```
