# Simulated Small Council pursuit data (bipartite)

A 10x7 matrix of simulated "pursuit intensity" scores representing how
aggressively each Great House from *Game of Thrones* seeks each Small
Council position. Generated with deliberate actor effects (some houses
are more ambitious overall), partner effects (some positions are more
sought after), and unique dyadic effects (house-specific affinities for
particular roles). Serves as a bipartite (two-mode) example where
senders and receivers come from different populations.

## Usage

``` r
data(small_council)
```

## Format

A 10x7 numeric matrix. Rows are Great Houses: Stark, Lannister,
Targaryen, Baratheon, Tyrell, Martell, Greyjoy, Arryn, Tully, Bolton.
Columns are Small Council positions: Hand, Coin, Whispers, Ships, War,
Law, Faith.

## Examples

``` r
data(small_council)
fit = srm(small_council)
summary(fit)
#> Social Relations Model - Variance Decomposition
#> ================================================== 
#> 
#> Component                Variance  % Total
#> ------------------------------------------ 
#> Actor                      2.0486    37.8%
#> Partner                    1.5662    28.9%
#> Unique                     1.7980    33.2%
```
