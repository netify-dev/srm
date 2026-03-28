# srm

## Overview

The `srm` package provides tools for analyzing network data via the Social Relations Model. It decomposes networks into actor (sender), partner (receiver), and unique dyadic components, computes variance partitions and covariance structures, and supports permutation-based inference. Works with both cross-sectional and longitudinal networks, including bipartite designs. Part of the [netify](https://github.com/netify-dev/netify) ecosystem.

## Installation

```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("netify-dev/srm")
```

## Usage

```r
library(srm)

# fit SRM to included classroom data
data(classroom)
fit = srm(classroom)
fit
summary(fit)

# visualize
plot(fit, type = "actor", n = 8)
plot(fit, type = "variance")

# permutation test
pt = permute_srm(fit, n_perms = 500, seed = 6886)
print(pt)
```

The package also works directly with `netify` objects:

```r
library(netify)
data(atop_EA)

net = netify::netify(
  input = atop_EA,
  actor1 = "country1",
  actor2 = "country2",
  symmetric = TRUE,
  sum_dyads = TRUE,
  diag_to_NA = TRUE,
  missing_to_zero = TRUE
)

fit = srm(net)
summary(fit)
```

See the vignettes for detailed examples: `srm-overview` (component-by-component walkthrough), `pipeline` (end-to-end analysis), `bipartite` (two-mode networks), and `methodology` (mathematical framework).

## Contributors

Cassy Dorff (Vanderbilt University), Shahryar Minhas (Michigan State University), Tosin Salau (Michigan State University)
