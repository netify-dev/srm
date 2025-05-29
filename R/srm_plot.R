# Define global variables used in ggplot
utils::globalVariables(c("effect", "id", "color", "actor", "partner"))

#' Create plots for SRM actor, dyadic, partner effects
#'
#' This function creates visualizations of social relations model effects including
#' actor effects, partner effects, and dyadic effects.
#'
#' @param mat A matrix or list of matrices from srm_effects()
#' @param type A string specifying the type of effect to plot: "actor", "partner", or "dyadic"
#' @param n The number of top actors/partners to plot
#' @param time A string specifying the time point to plot
#' @param facet A logical indicating whether to facet multiple matrices
#' @return A ggplot object or list of ggplot objects
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort 
#' @importFrom stats reorder
#' @examples
#' # Simple example with a basic matrix
#' test_matrix <- matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow=3, ncol=3)
#' rownames(test_matrix) <- colnames(test_matrix) <- c("A", "B", "C")
#' actor_effects <- srm_effects(test_matrix, type = "actor")
#' srm_plot(actor_effects, type = "actor", n = 3)
#' 
#' \dontrun{
#' # Example taken from GED data
#' # subsetted to Mexico with resulting df as a dyadic dataframe
#' 
#' UCDP_GED <- read_excel(paste0(pth,"GEDEvent_v24_1.xlsx"))
#' mexico <- UCDP_GED |>
#' filter(country=='Mexico')
#'
#' # turn dyadic df to a network object using netify using civ death as weight
#' # so we have a matrix for crossectional and list for longitudinal
#' 
#' # cross-sectional
#' 
#' mex_network_civ <- netify(
#' dyad_data = mexico,
#' actor1 = 'side_a',
#' actor2 = 'side_b',
#' weight = 'deaths_civilians',
#' symmetric = TRUE,
#' sum_dyads = TRUE,
#' diag_to_NA = TRUE,
#' missing_to_zero = TRUE
#' )
#' 
#' # longitudinal
#' mex_network_long <- netify(
#'   dyad_data = mexico,
#'   actor1 = 'side_a',
#'   actor2 = 'side_b',
#'   time = 'year',      
#'   weight = 'deaths_civilians',
#'   symmetric = TRUE,
#'   sum_dyads = TRUE,
#'   diag_to_NA = TRUE,
#'   missing_to_zero = TRUE
#'  )  
#' 
#' # plot actor effects
#' # before running, matrix should be from srm_effects
#' 
#' actor_effect <- srm_effects(mex_network_long, type = "actor")
#' 
#' actor_plot <- srm_plot(actor_effect, type = "actor", n = 5)
#' 
#' # plot dyadic effects
#' dyadic_effect <- srm_effects(mex_network_long, type = "unique")
#' dyadic_plot <- srm_plot(dyadic_effect, type = "dyadic", time= "2015")
#' 
#' # for longitudinal, can use facet argument for actor/partner effect
#' actor_plot_facet <- srm_plot(actor_effect, type = "actor", facet=TRUE)
#' }
#' @export
srm_plot <- function(mat, type = c("actor", "partner", "dyadic"), n = 10, time = NULL, facet = FALSE) {
  
  type <- match.arg(type)
  
  # Check if input is matrix or list of matrices
  if (!is.matrix(mat) && !is.list(mat)) {
    cli::cli_abort("Input must be a matrix or a list of matrices.\nYou can use the netify package to create a matrix from your data.")
  }
  
  # Convert matrix to list if needed
  if (is.matrix(mat)) {
    mat <- list(mat)
  }
  
  # Filter by time if requested
  if (!is.null(time)) {
    if (is.null(names(mat))) {
      cli::cli_abort("Time filtering requested, but list of matrices has no names.")
    }
    invalid_t <- time[!time %in% names(mat)]
    if (length(invalid_t) > 0) {
      cli::cli_abort(paste("Time point(s) not found:", paste(invalid_t, collapse = ", ")))
    }
    mat <- mat[time]
  }
  
  # Use facet for actor/partner effects with multiple matrices if requested
  if (facet && type %in% c("actor", "partner") && length(mat) > 1) {
    # combined data frame for all matrices
    combined_data <- NULL
    
    for (i in seq_along(mat)) {
      x <- mat[[i]]
      period_name <- names(mat)[i]
      if (is.null(period_name)) period_name <- paste("Matrix", i)
      
      # Extract data based on type
      if (type == "actor") {
        temp_data <- data.frame(
          id = rownames(x),
          effect = x[, 1],
          period = period_name
        )
      } else {
        temp_data <- data.frame(
          id = colnames(x),
          effect = x[1, ],
          period = period_name
        )
      }
      
      # Sort and limit to top n within each period
      temp_data <- temp_data[order(abs(temp_data$effect), decreasing = TRUE), ]
      temp_data <- temp_data[1:min(n, nrow(temp_data)), ]
      
      # Grayscale for effects
      temp_data$color <- ifelse(temp_data$effect > 0, "#252525", "#737373")
      
      # Combine with main df
      if (is.null(combined_data)) {
        combined_data <- temp_data
      } else {
        combined_data <- rbind(combined_data, temp_data)
      }
    }
    
    # Create the faceted plot
    p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = effect, y = reorder(id, effect), fill = color)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~ period, scales = "free", ncol = 2) +
      ggplot2::scale_fill_identity() +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 7),
        axis.text.x = ggplot2::element_text(size = 7, angle = 45, hjust = 1),
        plot.margin = ggplot2::unit(c(0.3, 0.5, 0.3, 0.5), "cm"),
        panel.spacing = ggplot2::unit(1, "cm"),
        strip.text = ggplot2::element_text(color = "white", size = 10),
        strip.background = ggplot2::element_rect(fill = "grey30", color = NA)
      ) +
      ggplot2::xlab("") +
      ggplot2::ylab("")
    
    return(p)
  }
  else if (facet && type == "dyadic" && length(mat) > 1) {
    
    cli::cli_warn("Faceting multiple dyadic heatmaps is not supported. Returning separate plots instead.")
  }
  
  # For single matrices, dyadic plots, or when faceting is false
  res <- lapply(seq_along(mat), function(i) {
    
    x <- mat[[i]]
    
    # Different plotting logic based on type
    if (type %in% c("actor", "partner")) {
      # Logic for actor/partner effects
      if (type == "actor") {
        ggdata <- data.frame(
          id = rownames(x),
          effect = x[, 1]
        )
      } else {
        ggdata <- data.frame(
          id = colnames(x),
          effect = x[1, ]
        )
      }
      
      # Sort by abs value
      ggdata <- ggdata[order(abs(ggdata$effect), decreasing = TRUE), ]
      
      # Top n actors
      n <- min(n, nrow(ggdata))
      ggdata <- ggdata[1:n, ]
      
      # Grayscale for effects
      ggdata$color <- ifelse(ggdata$effect > 0, "#252525", "#737373")
      
      # Create bar plot
      ggplot2::ggplot(ggdata, ggplot2::aes(x = effect, y = reorder(id, effect), fill = color)) +
        ggplot2::geom_col() +
        ggplot2::ggtitle(names(mat)[i]) +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = 7)
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab("")
    } 
    else if (type == "dyadic") {
      # First, filter the matrix to top n actors by total effect
      # Calculate total effect magnitude for each actor (sum of absolute row and column values)
      row_sums <- rowSums(abs(x))
      col_sums <- colSums(abs(x))
      total_effects <- row_sums + col_sums
      
      # Get top n actors
      top_n <- min(n, length(total_effects))
      top_actors <- names(sort(total_effects, decreasing = TRUE)[1:top_n])
      
      # Filter the matrix to only include top actors
      x_filtered <- x[top_actors, top_actors]
      
      # Reshape filtered data for ggplot
      ggdata <- reshape2::melt(x_filtered)
      colnames(ggdata) <- c("actor", "partner", "effect")
      
      # heatmap using geom_tile() with grayscale
      ggplot2::ggplot(ggdata, ggplot2::aes(x = actor, y = partner, fill = effect)) +
        ggplot2::geom_tile() +
        ggplot2::ggtitle(names(mat)[i]) +
        ggplot2::scale_fill_gradient2(low = "#f7f7f7", mid = "white", high = "#252525", midpoint = 0) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::xlab("") +
        ggplot2::ylab("")
    }
  })
  
  
  if (length(res) == 1) {
    return(res[[1]])
  } else {
    return(res)
  }
}