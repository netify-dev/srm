#' This function calculates actor, partner,and unique effects via srm framework
#'
#' @param mat A matrix or list of matrices
#' @param type A string specifying the type of effect to calculate.
#' @param time A string specifying the time point(s) to calculate the effect
#' @return A vector or list (if longitudinal) of actor/partner/unique effects
#' @importFrom netify netify 
#' 
#' 
#' @examples
#' # Simple example 
#' test_matrix <- matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow=3, ncol=3)
#' 
#' rownames(test_matrix) <- colnames(test_matrix) <- c("A", "B", "C")

#' actor_effect <- srm_effects(test_matrix, type = "actor")
#' 
#' \dontrun{
#' #Example taken from GED data
#' #subsetted to Mexico with resulting df as a dyadic dataframe
#' 
#' UCDP_GED <- read_excel(paste0(pth,"GEDEvent_v24_1.xlsx"))
#' mexico <- UCDP_GED |>
#' filter(country=='Mexico')

#' #turn dyadic df to a network object using netify using civ death as weight
#' #so we have a matrix for crossectional and list for longitudinal
#' 
#' #cross-sectional
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
#' #longitudinal

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
#' # Calculate actor effects
#' actor_effects <- srm_effects(mex_network_civ, type = "actor")
#'
#' # Calculate partner effects 
#' partner_effects <- srm_effects(mex_network_civ, type = "partner")
#'
#' # Calculate unique effects
#' unique_effects <- srm_effects(mex_network_civ, type = "unique")
#' 
#' #Calculate actor effects for a specific time point
#' actor_2015 <- srm_effects(mex_network_long, type = "actor", time = "2015")
#' 
#' #Calculate actor effects for  specific time points
#' yrs <- c("2015", "2016", "2017")
#' actor_long <- srm_effects(mex_network_long, type = "unique", time = yrs)
#' }
#' @export
srm_effects = function(mat, type= c("actor", "partner", "unique"), time = NULL ) {
  
  #check if input is a matrix or a list of matrices
  if (!is.matrix(mat) && !is.list(mat)) {
    cli::cli_abort("Input must be a matrix or a list of matrices.
                   You can use the netify package to create a matrix from your data.")
  }
  
  type <- match.arg(type)
  
  netify_warn <- FALSE
  
  if (is.matrix(mat)) {
    if (inherits(mat, "netify")) {
      cli::cli_alert_info("Converting 'netify' matrix to standard matrix.")
      class(mat) <- "matrix"
    }
    mat <- list(mat)
  } else if (is.list(mat)) {
    mat <- lapply(mat, function(m) {
      
      if (inherits(m, "netify")) {
        if (!netify_warn) {
          cli::cli_alert_info("Converting one or more 'netify' matrices in the list to standard matrices.")
          netify_warn <<- TRUE
        }
        class(m) <- "matrix"  # <- this happens for every netify matrix
      }
      return(m)
    })
  }
  
  #if matrix, turn into a list of matrices
  if (is.matrix(mat)) {
    mat <-list(mat)
  }
  
  # filter to requested time points if specified
  if (!is.null(time)) {
    if (is.null(names(mat))) {
      cli::cli_abort("Time filtering requested, but list of matrices has no names.")
    }
    # check time points exist
    invalid_t <- time[!time %in% names(mat)]
    if (length(invalid_t) > 0) {
      cli::cli_abort(paste("Time point(s) not found:", 
                           paste(invalid_t, collapse=", ")))
    }
    mat <- mat[time]
  }
  
  #srm function
  results = lapply(mat, function(x){
    if (dim(x)[1]!= dim(x)[2]) {
      cli::cli_abort("Matrix must be square")
    }
    
    # set diagonal to zero if not already
    diag(x) <- 0
    
    # extract dimensions
    n <- rd <- dim(x)[1]
    cd <- dim(x)[2]
    
    # calculate constants
    d1 <- (n - 1)  
    d2 <- (n - 2)  
    d3 <- (n * d2) 
    
    # calculate row, column, and total means
    x_r <- rowSums(x) / d1  
    x_c <- colSums(x) / d1  
    x_t <- sum(x) / (n * d1) 
    
    # actor effect (a_hat) represents how an individual i tends to rate others
    a_hat <- ((d1^2) / d3 * x_r) + ((d1 / d3) * x_c) - ((d1 / d2) * x_t)
    
    # Partner effect (b_hat) represents how an individual i tends to be perceived by others
    b_hat <- ((d1^2) / d3 * x_c) + ((d1 / d3) * x_r) - ((d1 / d2) * x_t)
    
    # calculate actor effect
    if (type == "actor"){
      
      a_hat_output <- as.matrix(a_hat, nrow=rd, ncol=1)
      rownames(a_hat_output) <- rownames(x)
      colnames(a_hat_output) <- list(c("actor effect for i"))
      
      # add class attribute
      class(a_hat_output) <- c("actor_effect", class(a_hat_output))
      
      return(a_hat_output)
      
    }
    
    # calculate partner effect
    else if (type == "partner"){
      
      
      b_hat_output <- as.matrix(b_hat, nrow=1, ncol=cd)
      rownames(b_hat_output) <- rownames(x)
      colnames(b_hat_output) <- list(c("partner effect for i"))
      
      # add class attribute
      class(b_hat_output) <- c("partner_effect", class(b_hat_output))
      
      return(b_hat_output)
      
    }
    
    # calculate unique effect
    else if (type == "unique"){
      
      # unique effect (g_hat) represents the unique effect of i on j
      
      g_hat<-matrix(NA, dim(x)[1], dim(x)[2])
      
      diag(g_hat)<-0.0  # enforce diagonal equal to zero
      
      diffMat= matrix(a_hat,nrow(x),ncol(x)) + 
        matrix(b_hat,nrow(x),ncol(x),byrow=T) + 
        matrix(x_t,nrow(x),ncol(x))
      
      g_hat= x-diffMat; diag(g_hat)=0
      
      g_hat_output<-as.matrix(g_hat, nrow=rd, ncol=cd) 
      rownames(g_hat_output)<-rownames(x)
      colnames(g_hat_output)<-rownames(x)
      
      # add class attribute
      class(g_hat_output) <- c("unique_effect", class(g_hat_output))
      
      return(g_hat_output)  
    }
    
  }
  )
  
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}


