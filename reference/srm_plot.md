# Create plots for SRM actor, dyadic, partner effects

Visualizes SRM effects computed by
[`srm_effects()`](https://netify-dev.github.io/srm/reference/srm_effects.md).
Actor and partner effects produce horizontal bar charts sorted by
absolute magnitude (dark bars = positive, gray bars = negative). Unique
(dyadic) effects produce a heatmap on a diverging color scale: blue
cells indicate weaker-than-expected ties, white cells are close to
predicted, and dark cells indicate stronger-than-expected ties.

## Usage

``` r
srm_plot(
  mat,
  type = c("actor", "partner", "dyadic"),
  n = 10,
  time = NULL,
  facet = FALSE
)
```

## Arguments

- mat:

  An effect matrix or list of effect matrices from
  [`srm_effects()`](https://netify-dev.github.io/srm/reference/srm_effects.md).
  Must not be an `srm` object.

- type:

  One of `"actor"`, `"partner"`, or `"dyadic"`. For bipartite networks,
  `"dyadic"` selects top rows and columns independently.

- n:

  Maximum number of actors/partners to display, selected by absolute
  magnitude (default 10).

- time:

  Optional character vector of time-point names to subset when `mat` is
  a named list.

- facet:

  If `TRUE` and `mat` is a list, facet actor/partner plots across time.
  Dyadic heatmaps are returned as separate plots instead.

## Value

A `ggplot` object, or a list of `ggplot` objects for unfaceted
longitudinal input.

## Details

For the unified interface (plot directly from an `srm` object), use
[`plot.srm()`](https://netify-dev.github.io/srm/reference/plot.srm.md)
instead.

## Examples

``` r
# actor effect bar plot
data(classroom)
actor_eff = srm_effects(classroom, type = "actor")
srm_plot(actor_eff, type = "actor", n = 8)


# dyadic heatmap
unique_eff = srm_effects(classroom, type = "unique")
srm_plot(unique_eff, type = "dyadic", n = 8)
```
