# Define global variables used in ggplot to avoid R CMD check notes
utils::globalVariables(c("effect", "id", "color", "actor", "partner"))

#' plot dyadic effect
#'
#' This function creates a heatmap of dyadic effects from srm_effects(unique)?. 
#' It colors negative values in red and positive values in blue 
#'
#' @param mat matrix from srm_effects(unique).
#' @return A ggplot heatmap object.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort 
#' @export
#' 
srm_plotdyadic<- function(mat) {
  
  # Check if input is a matrix.... force its from srm_effects(unique)?
  if (!is.matrix(mat)) {
    cli::cli_abort("Input must be a matrix.
                   You can use the netify package to create a matrix from your data")
  
  }
  
  # Reshape data for ggplot
  ggdata <- reshape2::melt(mat)
  colnames(ggdata) <- c("actor", "partner", "effect")
  
  # Create heatmap using `geom_tile()`
  ggplot(ggdata, aes(x = actor, y = partner, fill = effect)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}


#' plot actor effect
#' function plots actor and partner effect from srm_effects(actor/partner)
#' 
#' 
#' @param mat matrix from srm_effects(actor/partner).
#' @param n number of top actors/partners to plot.
#' @param type A string specifying the type of effect to plot.
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom cli cli_abort
#' @importFrom stats reorder 
#' @export
#' 

srm_plotactor= function(mat, n=10, type = c("actor", "partner")){
  
  type <- match.arg(type)
  
  # Check if input is a matrix
  if (!is.matrix(mat)) {
    cli::cli_abort("Input must be a matrix from srm_effects().")
  }
  
  #convert matrix to data frame
  if(type== "actor"){
    
    #actor effects
    
    ggdata = data.frame(
      id= rownames(mat),
      effect = mat[,1]
    )
  }
  
  else{
    #partner effects
    ggdata= data.frame(
      id =colnames(mat),
      effect = mat[1,]
    )
  }
  
  # Sort by abs value
  ggdata <- ggdata[order(abs(ggdata$effect), decreasing = TRUE), ]
  
  #top n actors
  n = min(n, nrow(ggdata))
  ggdata <- ggdata[1:n, ]
  
  #color
  ggdata$color <- ifelse(ggdata$effect > 0, "blue", "red")
  
  # Create plot
  
  ggplot(ggdata, aes(x = effect, y = reorder(id, effect), fill = color)) +
    geom_col() +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 7)
    )
}






