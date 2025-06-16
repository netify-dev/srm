# srm

## Overview <img src="https://github.com/netify-dev/srm/blob/main/man/figures/srm_hex.png" align = "right"  alt = "hex" width="200px"> 

The proposed `srm` package proposes tools for analyzing nodal dyadic relationships via the social relations model. It breaks down and analyzes networks into nodal and dyadic dependencies and allows users to understand actor, partner, and dyadic effects in network data.

## Installation

    if (!require(devtools)) {
        install.packages("devtools")
      }
      library(devtools)

      install_github("netify-dev/srm")
      
## Usage

See our `srm-overview` vignette for detailed examples and documentation. To get started, supply network data (matrix or netify object) to the `srm` functions. You can analyze actor effects, partner effects, unique dyadic relationships, and generate summary statistics. For example, to analyze actor effects, we use the code below:


    library(srm)
    data("atop_EA")
    
    # transform data into matrix or netify object
    dat <- netify::netify(
      dyad_data = atop_EA,
      actor1 = 'country1',
      actor2 = 'country2',
      symmetric = TRUE,
      sum_dyads = TRUE,
      diag_to_NA = TRUE,
      missing_to_zero = TRUE
    )
    
    # calculate actor effects to identify actor patterns
    actor_effects <- srm_effects(dat, type = "actor")
    
    # generate summary statistics
    row_means <- srm_stats(dat, type = "rowmeans")
    
    # visualize the results
    srm_plot(actor_effects, type = "actor", n = 8)


## Contributors 

Cassy Dorff (Vanderbilt University), Shahryar Minhas (Michigan State University), Tosin Salau (Michigan State University)
