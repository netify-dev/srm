# Permutation test for SRM variance components

Tests whether actor, partner, and unique variance components are
significantly different from zero using a permutation (randomization)
approach that preserves network size but destroys relational structure.

## Usage

``` r
permute_srm(mat, n_perms = 1000, time = NULL, seed = NULL)
```

## Arguments

- mat:

  A square matrix, list of square matrices, `netify` object, or an
  object of class `srm`. Bipartite (rectangular) networks are not
  supported; an error is raised if the input is non-square.

- n_perms:

  Number of permutations (default 1000).

- time:

  Optional time points to subset to.

- seed:

  Optional random seed for reproducibility.

## Value

An object of class `srm_permtest` containing:

- observed:

  Named vector of observed variance components.

- perm_dist:

  Matrix of permuted variance components (`n_perms` x components).

- p_values:

  One-tailed p-values for each component.

- n_perms:

  Number of permutations used.

## Details

For each permutation, the rows and columns of the adjacency matrix are
independently permuted, destroying actor and partner structure while
preserving the marginal distribution. The SRM variance components are
recomputed on each permuted matrix and the observed values are compared
to this null distribution.

**Interpreting unique variance p-values:** Permutation typically
increases dyadic noise, so the observed unique variance is usually
smaller than the null distribution. A p-value of 1.0 for unique variance
is expected and does not imply that dyadic effects are absent.

## Examples

``` r
# \donttest{
data(classroom)
pt = permute_srm(classroom, n_perms = 200, seed = 6886)
print(pt)
#> SRM Permutation Test
#> Permutations: 200
#> -------------------------------------------------- 
#> Component                Observed Mean(Null)      p
#> -------------------------------------------------- 
#> Actor Var                  1.4000     1.1262  0.045 *
#> Partner Var                0.8539     0.6584  0.045 *
#> Unique Var                 1.0546     2.2181  1.000 
#> Relationship Cov           0.6885     0.0062  0.020 *
#> Actor-Partner Cov         -0.2525     0.0366  0.806 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1
# }
```
