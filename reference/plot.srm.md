# Plot an srm object

Dispatches to actor, partner, dyadic, or variance partition plots.

## Usage

``` r
# S3 method for class 'srm'
plot(
  x,
  type = c("actor", "partner", "dyadic", "variance"),
  n = 10,
  time = NULL,
  facet = FALSE,
  ...
)
```

## Arguments

- x:

  An object of class `srm`.

- type:

  One of `"actor"`, `"partner"`, `"dyadic"`, or `"variance"`.

- n:

  Maximum number of actors to display (for actor/partner/dyadic).

- time:

  Optional time point(s) to restrict to.

- facet:

  Logical; facet across time for actor/partner plots.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object or list of `ggplot` objects.

## Examples

``` r
data(classroom)
fit = srm(classroom)
plot(fit, type = "actor")

plot(fit, type = "variance")

```
