#' Density of the Dirichlet-Multinomial Distribution
#'
#' @description Density function for the Dirichlet-multinomial distribution.
#' @details Computes the probability mass of the Dirichlet-multinomial distribution. Under the proportion parameterization, the alpha parameters of the conventional Dirichlet-multinomial distribution are derived as the product of a proportion vector (`p`) and an exponentiated precision parameter ([`exp`][exp()]`(theta)`). The precision parameter controls the degree of overdispersion relative to the multinomial distribution, where higher values of `theta` are associated with reduced overdispersion. When `theta = 0`, the alpha parameters of the conventional Dirichlet-multinomial distribution are equal to the proportion vector (`p`). To ensure a simplex, the values of `p` (vector or matrix records) are internally normalized to sum to one. If `alpha` is provided, then the conventional alpha parameterization of the Dirichlet-multinomial distribution is used.
#' @param x Numeric vector or matrix of counts. If matrix, then a vector of probability densities is returned with an element for each record of the matrix. Matrix records represent observations, and matrix fields represent dimensions.
#' @param p Numeric vector or matrix of proportions. Matrix records represent observations, and matrix fields represent dimensions. If vector, then `p` is recycled for each record of matrix `x`.
#' @param theta Numeric vector. Precision parameter with domain (`-Inf`, `Inf`). If scalar, then `theta` is recycled for each record of matrix `x`.
#' @param alpha Numeric vector or matrix of conventional alpha values. Matrix records represent observations, and matrix fields represent dimensions. If vector, then `alpha` is recycled for each record of matrix `x`.
#' @param log Logical scalar. If `TRUE`, then probabilities are given as `log(density)`.
#' @returns Numeric vector of probability densities.
#' @seealso
#' [`waic`][waic()] for generic function to compute widely applicable information criterion. \cr \cr
#' [`dmWAIC`][dmWAIC()] for computing widely applicable information criteria for Dirichlet-multinomial regression models.
#' @examples
#' # Compute log probability density.
#' ddirmult(x=c(33,115,95,359),
#'          p=c(0.075,0.201,0.175,0.549),
#'          theta=4.027,log=TRUE)
#' @export
ddirmult<-function(x,p,theta,alpha,log=FALSE){
  
  # Check x parameter.
  ## If x is a vector.
  if(is.vector(x)){
    ## Convert vector to matrix.
    x<-matrix(x,nrow=1,byrow=FALSE)
  }
  ## Check that x is a matrix.
  if(!is.matrix(x)) stop("x must be a vector or matrix.")
  ## Check that x is numeric.
  if(!is.numeric(x)) stop("x must be numeric integer.")
  ## Check that x is integer.
  if(any(x %% 1 != 0)) stop("x must be numeric integer.")
  ## Check that x is positive.
  if(any(x < 0)) stop("x must be positive.")
  
  # If alpha is missing.
  if(missing(alpha)){
    
    # Check p parameter.
    ## If p is a vector.
    if(is.vector(p)){
      ## Convert vector to matrix.
      p<-matrix(p,nrow=1,byrow=FALSE)
      ## Match dimensions with x.
      p<-p[rep(x=1,times=nrow(x)),,drop=FALSE]
    }
    ## Check that p is a matrix.
    if(!is.matrix(p)) stop("p must be a vector or matrix.")
    ## Check that p has proper dimensions.
    if(!identical(dim(p),dim(x))) stop("p has improper dimensions.")
    ## Check that p is numeric.
    if(!is.numeric(p)) stop("p must be numeric.")
    ## Check that p is positive.
    if(any(p < 0)) stop("p must be positive.")
    
    # Check theta parameter.
    ## Check that theta is a vector.
    if(!is.vector(theta)) stop("theta must be a vector.")
    ## If theta has length one.
    if(length(theta)==1){
      ## Match length with x.
      theta<-rep(x=theta,times=nrow(x))
    }
    ## Check that theta has proper length.
    if(length(theta)!=nrow(x)) stop("theta has improper length.")
    ### Check that theta is numeric.
    if(!is.numeric(theta)) stop("theta must be numeric.")
    
    # Normalize p to ensure simplex.
    p<-normalize(x=p)
    
    # Derive alpha parameters.
    alpha<-p*exp(theta)
    
  }else{ # If alpha is provided.
    
    # Check alpha parameter.
    ## If alpha is a vector.
    if(is.vector(alpha)){
      ## Convert vector to matrix.
      alpha<-matrix(alpha,nrow=1,byrow=FALSE)
      ## Match dimensions with x.
      alpha<-alpha[rep(x=1,times=nrow(x)),,drop=FALSE]
    }
    ## Check that alpha is a matrix.
    if(!is.matrix(alpha)) stop("alpha must be a vector or matrix.")
    ## Check that alpha has proper dimensions.
    if(!identical(dim(alpha),dim(x))) stop("alpha has improper dimensions.")
    ## Check that alpha is numeric.
    if(!is.numeric(alpha)) stop("alpha must be numeric.")
    ## Check that alpha is positive.
    if(any(alpha < 0)) stop("alpha must be positive.")
    
    # Check for other parameters.
    ## Check if p parameter is provided.
    if(!missing(p)) warning("p is ignored when alpha is provided.")
    ## Check if theta parameter is provided.
    if(!missing(theta)) warning("theta is ignored when alpha is provided.")
    
  }
  
  # Check log parameter.
  ## Check that log is a vector.
  if(!is.vector(log)) stop("log must be a vector.")
  ## Check that log is logical.
  if(!is.logical(log)) stop("log must be logical.")
  ## Check that log has length 1.
  if(length(log)!=1) stop("log must have length 1.")
  ## Check that log is not NA.
  if(is.na(log)) stop("log cannot be NA.")
  
  # Define log probability mass function.
  lpmf<-function(x,alpha){
    
    # Calculate log probability mass.
    m<-lgamma(sum(alpha))+lgamma(sum(x)+1)-lgamma(sum(x)+sum(alpha))+
      sum(lgamma(x+alpha))-sum(lgamma(alpha))-sum(lgamma(x+1))
    
    # Return log probability mass.
    return(m)
    
  }
  
  # Vectorize log probability calculations.
  d<-sapply(X=1:nrow(x),FUN=function(i) lpmf(x=x[i,],alpha=alpha[i,]))
  
  # Exponentiate if log is FALSE.
  if(!log) d<-exp(d)
  
  # Return density.
  return(d)
  
}