#' This function calculates srm statistics
#' 
#' @param mat A matrix or list of matrices.
#' @param type A string specifying the type of statistic to calculate.
#' @param time A string specifying the time point(s) to calculate.
#' @return A vector or list (if longitudinal) of calculated network stats 
#' @importFrom cli cli_abort
#' @examples
#' # Simple example 
#' test_matrix <- matrix(c(0, 1, 2, 2, 0, 1, 1, 2, 0), nrow=3, ncol=3)
#' 
#' rownames(test_matrix) <- colnames(test_matrix) <- c("A", "B", "C")

#' actor_totalmeans <- srm_stats(test_matrix, type = "totalmeans")
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
#' # Calculate srm stats for rowmeans
#' actor_mean <- srm_stats(mex_network_civ, type = "rowmeans")
#' 
#' # Calculate srm stats for a specific time point
#' 
#' actor_mean <- srm_stats(mex_network_long, type = "partner_var", time = "2015")
#' 
#' }
#' @export

srm_stats <-function(mat, type =c("rowmeans",
                                  "colmeans",
                                  "totalmeans",
                                  "actor_var",
                                  "unique_var", 
                                  "partner_var",
                                  "relationship_cov", 
                                  "actor_partner_cov"), time =NULL){
  
  # check if input is matrix or list of matrices
  if (!is.matrix(mat) && !is.list(mat)) {
    cli::cli_abort("Input must be a matrix or a list of matrices.\nYou can use the netify package to create a matrix from your data.")
  }
  
  type <- match.arg(type)
  
  # if matrix, wrap in list
  if (is.matrix(mat)) {
    mat <- list(mat)
  }
  
  # filter by time (if specified)
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
  
  
  #lapply to calculate the stats
  results<-lapply(mat, function(x){
    if (dim(x)[1] != dim(x)[2]) {
      cli::cli_abort("Matrix must be square")
    }
    
    
    # Set diagonal to zero if not already
    diag(x) <- 0
    
    # Extract dimensions
    n <- rd <- dim(x)[1]
    cd <- dim(x)[2]
    dim<-(dim(x))
    
    # Calculate constants
    d1 <- (n - 1)  
    d2 <- (n - 2)  
    d3 <- (n * d2) 
    d4 <- ((d1) * (d2))
    d5 <- (((d1) * (d2) / (2)) - 1)
    
    #calculate row, column, and total means
    x_r <- rowSums(x) / d1
    x_c <- colSums(x) / d1
    x_t <- sum(x) / (n * d1)
    
    if (type == "rowmeans") {
      names(x_r) <- rownames(x)
      return(x_r)
    }
    else if (type == "colmeans") {
      names(x_c) <- colnames(x)
      return(x_c)
    }
    else if (type == "totalmeans") {
      return(x_t)
    }
    
    a_hat<-((d1^2)/(d3) * x_r) + (((d1)/(d3))*x_c) - ((d1)/(d2) *x_t) #actor i perceives others
    b_hat<-((d1^2)/(d3) * x_c) + (((d1)/(d3))*x_r) - ((d1)/(d2) * x_t) #actor i is perceieved 
    
    g_hat<-matrix(NA, dim(x)[1], dim(x)[2])
    diag(g_hat)<-0.0  # enforce diagonal equal to zero
    
    
    diffMat=matrix(a_hat,nrow(x),ncol(x))+matrix(b_hat,nrow(x),ncol(x),byrow=T)+matrix(x_t,nrow(x),ncol(x))
    g_hat=x-diffMat; diag(g_hat)=0
    
    g_hatdiffs<-(g_hat - t(g_hat))^2 #all the differences
    diag(g_hatdiffs)<-0.0
    g_hatsums<-((g_hat + t(g_hat))/2)^2 #dividing by 2 is in equation #a_hat sums twice as big as g_ij + h_ji /2 
    diag(g_hatsums)<-0.0
    
    q1<-sum(g_hatsums) #summing squares of g_hat sum
    q2<-sum(g_hatdiffs)/2
    
    s_2g<-((q1/d5) + (q2/d4))/2 #took out *1/2
    s_gg<-((q1/d5) - (q2/d4))/2
    
    if(type== "actor_var"){
      
      s_2a<-((sum(a_hat^2)/(d1)) - ((s_2g*(d1))/(d3)) - ((s_gg)/(d3)))
      
      return(s_2a)
    }
    
    else if(type== "unique_var"){
      
      return(s_2g)
    }
    
    else if(type== "partner_var"){
      
      s_2b<-((sum(b_hat^2)/(d1)) - ((s_2g*(d1))/(d3)) - ((s_gg)/(d3)))
      
      return(s_2b)
    }
    
    else if(type== "relationship_cov"){
      
      return(s_gg)
    }
    
    else if(type== "actor_partner_cov"){
      
      s_ab<-((sum(a_hat * b_hat))/(d1)) - ((s_gg*(d1))/(d3)) - (s_2g/(d3))
      
      return(s_ab)
    }
    
    
  })
  
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

