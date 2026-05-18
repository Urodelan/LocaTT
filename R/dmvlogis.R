#' Density of the Multivariate Logistic Distribution
#'
#' @description Density function for the multivariate logistic distribution.
#' @details Computes the probability density of the multivariate logistic distribution. The multivariate logistic distribution is constructed using a Gaussian copula with logistic marginals. The probability density is the product of the densities of the logistic marginals, which is further multiplied by the density of a Gaussian copula of the transformed standard uniform margins (*i.e.*, probability integral transformation of the logistic marginals with [`stats::plogis`][stats::plogis()]).
#' @param x Numeric vector or matrix. Values of logistically-distributed marginals. If matrix, then a vector of probability densities is returned with an element for each record of the matrix. Matrix records represent observations, and matrix fields represent dimensions.
#' @param location Numeric vector. Location parameters of the logistic distribution. If `x` is a matrix, then `location` is recycled for each record of matrix `x`.
#' @param scale Numeric vector. Scale parameters of the logistic distribution. If `x` is a matrix, then `scale` is recycled for each record of matrix `x`.
#' @param R Numeric correlation matrix. If `x` is a matrix, then `R` is recycled for each record of matrix `x`.
#' @param log Logical scalar. If `TRUE`, then probabilities are given as `log(density)`.
#' @returns Numeric vector of probability densities.
#' @seealso
#' [`stats::dlogis`][stats::dlogis()] for density of the logistic distribution. \cr \cr
#' [`dcopula`][dcopula()] for density of the Gaussian copula.
#' @references Decani JS, and Stine RA. 1986. A note on deriving the information matrix for a logistic distribution. *The American Statistician*, 40(3): 220-222. DOI: 10.2307/2684541  \cr \cr
#' Song P. 2000. Multivariate dispersion models generated from Gaussian copula. *Scandinavian Journal of Statistics*, 27(2): 305-320. DOI: 10.1111/1467-9469.00191
#' @examples
#' # Define logistic margins.
#' x<-c(0.055,-1.625,0.329,-5.765)
#' 
#' # Define location parameters.
#' location<-c(0.477,-0.998,-0.776,0.064)
#' 
#' # Define scale parameters.
#' scale<-c(0.574,1.314,0.460,1.393)
#' 
#' # Define correlation matrix.
#' R<-matrix(data=c(1.000,-0.80,0.64,-0.512,
#'                  -0.800,1.00,-0.80,0.640,
#'                  0.640,-0.80,1.00,-0.800,
#'                  -0.512,0.64,-0.80,1.000),
#'            ncol=4,byrow=TRUE)
#' 
#' # Compute log probability density.
#' dmvlogis(x=x,location=location,
#'          scale=scale,R=R,
#'          log=TRUE)
#' @export
dmvlogis<-function(x,location,scale,R,log=FALSE){
  
  # Check x parameter.
  ## If x is a vector.
  if(is.vector(x)){
    ## Convert vector to matrix.
    x<-matrix(x,nrow=1,byrow=FALSE)
  }
  ## Check that x is a matrix.
  if(!is.matrix(x)) stop("x must be a vector or matrix.")
  ## Check that x is numeric.
  if(!is.numeric(x)) stop("x must be numeric.")
  
  # Check location parameter.
  ## Check that location is a vector.
  if(!is.vector(location)) stop("location must be a vector.")
  ## Check that location is numeric.
  if(!is.numeric(location)) stop("location must be numeric.")
  ## Check that location has proper dimensions.
  if(length(location)!=ncol(x)) stop("location has improper dimensions.")
  
  # Check scale parameter.
  ## Check that scale is a vector.
  if(!is.vector(scale)) stop("scale must be a vector.")
  ## Check that scale is numeric.
  if(!is.numeric(scale)) stop("scale must be numeric.")
  ## Check that scale is positive.
  if(any(scale < 0,na.rm=TRUE)) stop("scale must be positive.")
  ## Check that scale has proper dimensions.
  if(length(scale)!=ncol(x)) stop("scale has improper dimensions.")
  
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
  ## Check that R has proper dimensions.
  if(ncol(R)!=ncol(x)) stop("R has improper dimensions.")
  
  # Check log parameter.
  ## Check that log is a vector.
  if(!is.vector(log)) stop("log must be a vector.")
  ## Check that log is logical.
  if(!is.logical(log)) stop("log must be logical.")
  ## Check that log has length 1.
  if(length(log)!=1) stop("log must have length 1.")
  ## Check that log is not NA.
  if(is.na(log)) stop("log cannot be NA.")
  
  # Compute log probability density of logistic marginals.
  d<-colSums(apply(X=x,MARGIN=1,FUN=stats::dlogis,
                   location=location,scale=scale,
                   log=TRUE))
  
  # Convert to quantiles of the logistic distribution.
  u<-t(apply(X=x,MARGIN=1,FUN=stats::plogis,location=location,scale=scale))
  
  # Increment log probability density by Gaussian copula.
  d<-d+dcopula(u=u,R=R,log=TRUE)
  
  # Exponentiate if log is FALSE.
  if(!log) d<-exp(d)
  
  # Return density.
  return(d)
  
}