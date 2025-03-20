#' This function calculates actor, partner,and unique effects via srm framework
#'
#' @param x can be a df or matrix. if df- converted to matrix by netify
#' @param type A character string specifying the type of effect to calculate.
#' @return A vector of actor/partner/unique effects
#' @importFrom netify netify 
#' @export
srm_effects <- function(x, type= c("actor", "partner", "unique")) {
  
  #sociomatrix from netify
  
  if (!is.matrix(x)) {
    cli::cli_abort("Input must be a matrix.
                   You can use the netify package to create a matrix from your data")
    
  }
  
  type <- match.arg(type)
  
  # Check if matrix is square
  if (dim(x)[1] != dim(x)[2]) {
    stop("Matrix must be square")
  }
  
  # Set diagonal to zero if not already
  diag(x) <- 0
  
  # Extract dimensions
  n <- rd <- dim(x)[1]
  cd <- dim(x)[2]
  
  # Calculate constants
  d1 <- (n - 1)  
  d2 <- (n - 2)  
  d3 <- (n * d2) 
  
  # Calculate row, column, and total means
  x_r <- rowSums(x) / d1  
  x_c <- colSums(x) / d1  
  x_t <- sum(x) / (n * d1) 
  
  # The actor effect (a_hat) represents how an individual i tends to rate others
  a_hat <- ((d1^2) / d3 * x_r) + ((d1 / d3) * x_c) - ((d1 / d2) * x_t)
  
  # Partner effect (b_hat) represents how an individual i tends to be perceived by others
  b_hat <- ((d1^2) / d3 * x_c) + ((d1 / d3) * x_r) - ((d1 / d2) * x_t)
  
  # Calculate actor effect
  if (type == "actor"){
    
    a_hat_output <- as.matrix(a_hat, nrow=rd, ncol=1)
    rownames(a_hat_output) <- rownames(x)
    colnames(a_hat_output) <- list(c("actor effect for i"))
    
    return(a_hat_output)
    
  }
  
  # Calculate partner effect
  else if (type == "partner"){

  b_hat_output <- as.matrix(b_hat, nrow=1, ncol=cd)  
  rownames(b_hat_output) <- rownames(x)
  colnames(b_hat_output) <- list(c("partner effect for i"))
  
  return(b_hat_output)
  
  }
  
  # Calculate unique effect
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
  
    return(g_hat_output)  
  }
  
}



#' function for means, variance
#' 
#' 
#' @param x A square matrix 
#' @param type A character string specifying the type of effect to calculate.
#' @return A vector of actor/partner/unique effects
#' @export

srm_stats <-function(x, type =c("rowmeans",
                                "colmeans",
                                "totalmeans",
                                "actor_variance",
                                "dyadic/unique_variance", 
                                "partner_variance",
                                "dyadic/relationship_covariance", 
                                "actor_partner_covariance")){
  
  #sociomatrix from netify 
  
  type <- match.arg(type)
  
  # Check if matrix is square
  if (dim(x)[1] != dim(x)[2]) {
    stop("Matrix must be square")
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
  
  if(type== "actor_variance"){
    
    s_2a<-((sum(a_hat^2)/(d1)) - ((s_2g*(d1))/(d3)) - ((s_gg)/(d3)))
    
    return(s_2a)
  }
  
  else if(type== "dyadic/unique_variance"){
    
    return(s_2g)
  }
  
  else if(type== "partner_variance"){
    
    s_2b<-((sum(b_hat^2)/(d1)) - ((s_2g*(d1))/(d3)) - ((s_gg)/(d3)))
    
    return(s_2b)
  }
  
  else if(type== "dyadic/relationship_covariance"){
    
    return(s_gg)
  }
  
  else if(type== "actor_partner_covariance"){
    
    s_ab<-((sum(a_hat * b_hat))/(d1)) - ((s_gg*(d1))/(d3)) - (s_2g/(d3))
    
    return(s_ab)
  }
}




# 
#' This function calculates actor, partner, and unique effects via srm framework
#'
#' @param x A matrix or list of matrices (for longitudinal data)
#' @param type A character string specifying the type of effect to calculate.
#' @param time For longitudinal data (list of matrices), specify which time points to analyze
#' @return A matrix or list of matrices with actor/partner/unique effects
#' @importFrom netify netify
#' @export

srm_effects_long<- function(x, type = c("actor", "partner", "unique"), time = NULL) {
  
  # Match argument for type parameter
  type <- match.arg(type)
  
# Handle list input (for longitudinal data) has to be a list and matrix
  if (is.list(x) && !is.matrix(x)) {
  
  # Filter to requested time points if specified
    if (!is.null(time)) {
    # Check time points exist
    invalid_t <- time[!time %in% names(x)]
      if (length(invalid_t) > 0) {
      cli::cli_abort(paste("Time point(s) not found:", 
                           paste(invalid_t, collapse=", ")))
    }
    x <- x[time]
  }
  
  # Process each time point
  results <- list()
  for (t in names(x)) {
    results[[t]] <- srm_effects_long(x[[t]], type)
  }
  return(results)
}

  if (!is.matrix(x)) {
    cli::cli_abort("Input must be a matrix or a list of matrices.
                   You can use the netify package to create a matrix from your data")
  }
  
  # Check if matrix is square
  if (dim(x)[1] != dim(x)[2]) {
    stop("Matrix must be square")
  }
  
  # Set diagonal to zero if not already
  diag(x) <- 0
  
  # Extract dimensions
  n <- rd <- dim(x)[1]
  cd <- dim(x)[2]
  
  # Calculate constants
  d1 <- (n - 1)  
  d2 <- (n - 2)  
  d3 <- (n * d2) 
  
  # Calculate row, column, and total means
  x_r <- rowSums(x) / d1  
  x_c <- colSums(x) / d1  
  x_t <- sum(x) / (n * d1) 
  
  # The actor effect (a_hat) represents how an individual i tends to rate others
  a_hat <- ((d1^2) / d3 * x_r) + ((d1 / d3) * x_c) - ((d1 / d2) * x_t)
  
  # Partner effect (b_hat) represents how an individual i tends to be perceived by others
  b_hat <- ((d1^2) / d3 * x_c) + ((d1 / d3) * x_r) - ((d1 / d2) * x_t)
  
  # Calculate actor effect
  if (type == "actor"){
    
    a_hat_output <- as.matrix(a_hat, nrow=rd, ncol=1)
    rownames(a_hat_output) <- rownames(x)
    colnames(a_hat_output) <- list(c("actor effect for i"))
    
    return(a_hat_output)
    
  }
  
  # Calculate partner effect
  else if (type == "partner"){
    
    b_hat_output <- as.matrix(b_hat, nrow=1, ncol=cd)  
    rownames(b_hat_output) <- rownames(x)
    colnames(b_hat_output) <- list(c("partner effect for i"))
    
    return(b_hat_output)
    
  }
  
  # Calculate unique effect
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
    
    return(g_hat_output)  
  }
  
  
}






