#' Density of the Gaussian Copula
#'
#' @description Density function for the Gaussian copula.
#' @details Computes the probability density of the Gaussian copula. Given uniformly-distributed margins `u` on the interval \[`0`, `1`\], applies the inverse cumulative distribution function of the standard normal (*i.e.*, [`stats::qnorm`][stats::qnorm()]) to map uniform margins to normal scores. Then, uses equation 1 of Song (2000) with the normal scores to calculate the probability density of the Gaussian copula.
#' @param u Numeric vector or matrix of uniformly-distributed margins on the interval \[`0`, `1`\]. If matrix, then a vector of probability densities is returned with an element for each record of the matrix. Matrix records represent observations, and matrix fields represent dimensions.
#' @param R Numeric correlation matrix. If `u` is a matrix, then `R` is recycled for each record of matrix `u`.
#' @param log Logical scalar. If `TRUE`, then probabilities are given as `log(density)`.
#' @returns Numeric vector of probability densities.
#' @seealso
#' [`dmvlogis`][dmvlogis()] for density of the multivariate logistic distribution.
#' @references Song P. 2000. Multivariate dispersion models generated from Gaussian copula. *Scandinavian Journal of Statistics*, 27(2): 305-320. DOI: 10.1111/1467-9469.00191
#' @examples
#' # Define uniform margins.
#' u<-c(0.324,0.383,0.917,0.015)
#' 
#' # Define correlation matrix.
#' R<-matrix(data=c(1.000,-0.80,0.64,-0.512,
#'                  -0.800,1.00,-0.80,0.640,
#'                  0.640,-0.80,1.00,-0.800,
#'                  -0.512,0.64,-0.80,1.000),
#'            ncol=4,byrow=TRUE)
#' 
#' # Compute log probability density.
#' dcopula(u=u,R=R,log=TRUE)
#' @export
dcopula<-function(u,R,log=FALSE){
  
  # Check u parameter.
  ## If u is a vector.
  if(is.vector(u)){
    ## Convert vector to matrix.
    u<-matrix(u,nrow=1,byrow=FALSE)
  }
  ## Check that u is a matrix.
  if(!is.matrix(u)) stop("u must be a vector or matrix.")
  ## Check that u is numeric.
  if(!is.numeric(u)) stop("u must be numeric.")
  ## Check that u is between 0 and 1.
  if(any((u < 0) | (u > 1),na.rm=TRUE)) stop("u must be between 0 and 1.")
  
  # Check R parameter.
  ## Check that R is a matrix.
  if(!is.matrix(R)) stop("R must be a matrix.")
  ## Check that R is numeric.
  if(!is.numeric(R)) stop("R must be numeric.")
  ## Check that all R diagonals equal 1.
  if(any(diag(R)!=1)) stop("R diagonals must be 1.")
  ## Check that R is symmetric.
  if(!isSymmetric(R)) stop("R must be symmetric.")
  ## Check that R is between -1 and 1.
  if(any(abs(R) > 1,na.rm=TRUE)) stop("R must be between -1 and 1.")
  
  # Check that R has proper dimensions.
  if(ncol(R)!=ncol(u)) stop("R has improper dimensions.")
  
  # Check log parameter.
  ## Check that log is a vector.
  if(!is.vector(log)) stop("log must be a vector.")
  ## Check that log is logical.
  if(!is.logical(log)) stop("log must be logical.")
  ## Check that log has length 1.
  if(length(log)!=1) stop("log must have length 1.")
  ## Check that log is not NA.
  if(is.na(log)) stop("log cannot be NA.")
  
  # Map to standard normal space.
  z<-stats::qnorm(p=u,mean=0,sd=1)
  
  # Define log probability density function.
  lpdf<-function(z,R){
    
    # Calculate log probability density.
    d<- -0.5 * (log(det(R)) + t(z) %*% (solve(R)-diag(ncol(R))) %*% z)
    
    # Return log probability density.
    return(d)
    
  }
  
  # Vectorize log probability calculations.
  d<-apply(X=z,MARGIN=1,FUN=lpdf,R=R)
  
  # Exponentiate if log is FALSE.
  if(!log) d<-exp(d)
  
  # Return density.
  return(d)
  
}