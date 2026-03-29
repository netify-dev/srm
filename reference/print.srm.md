# Print an srm object

Print an srm object

## Usage

``` r
# S3 method for class 'srm'
print(x, ...)
```

## Arguments

- x:

  An object of class `srm`.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns `x`.

## Examples

``` r
mat = matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), 3, 3)
rownames(mat) = colnames(mat) = c("A", "B", "C")
fit = srm(mat)
print(fit)
#> Social Relations Model
#> ---------------------------------------- 
#> Mode:       unipartite
#> Actors:     3
#> Grand mean: 1.5000
#> 
#> Variance Decomposition:
#>   (Not estimable -- network too small, need n >= 4)
```
