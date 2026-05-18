#' Density of a Joint Species Distribution Model
#'
#' @description Density function for a joint species distribution model.
#' @details Computes the probability density of a joint species distribution model. The probability of observing a community is calculated as the product of the probabilities of observing each species. Observations for each species are Bernoulli-distributed, and species-specific probability densities are computed with [`stats::dbinom`][stats::dbinom()].
#' @param x Numeric vector or matrix. Binary values of species occurrence. If matrix, then a vector of probability densities is returned with an element for each record of the matrix. Matrix records represent sites, and matrix fields represent species.
#' @param psi Numeric vector or matrix. Probabilities of site occupancy. Matrix records represent sites, and matrix fields represent species. If vector, then `psi` is recycled for each record of matrix `x`.
#' @param log Logical scalar. If `TRUE`, then probabilities are given as `log(density)`.
#' @returns Numeric vector of probability densities.
#' @seealso
#' [`stats::dbinom`][stats::dbinom()] for density of the binomial distribution. \cr \cr
#' [`mlWAIC`][mlWAIC()] for computing widely applicable information criteria for joint species distribution models.
#' @references Wilkinson DP, Golding N, Guillera‐Arroita G, Tingley R, and McCarthy MA. 2021. Defining and evaluating predictions of joint species distribution models. *Methods in Ecology and Evolution*, 12(3): 394-404. DOI: 10.1111/2041-210X.13518
#' @examples
#' # Define species occurrence.
#' x<-c(1,0,0,1)
#' 
#' # Define occupancy probabilities.
#' psi<-c(0.886,0.391,0.139,0.991)
#' 
#' # Compute log probability density.
#' djsdm(x=x,psi=psi,log=TRUE)
#' @export
djsdm<-function(x,psi,log=FALSE){
  
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
  ## Check that x is binary.
  if(!all((x %in% 0:1) | is.na(x))) stop("x must be binary.")
  
  # Check psi parameter.
  ## If psi is a vector.
  if(is.vector(psi)){
    ## Convert vector to matrix.
    psi<-matrix(psi,nrow=1,byrow=FALSE)
    ## Match length with x.
    psi<-psi[rep(x=1,times=nrow(x)),,drop=FALSE]
  }
  ## Check that psi is a matrix.
  if(!is.matrix(psi)) stop("psi must be a vector or matrix.")
  ## Check that psi is numeric.
  if(!is.numeric(psi)) stop("psi must be numeric.")
  ## Check that psi is between 0 and 1.
  if(any((psi < 0) | (psi > 1),na.rm=TRUE)) stop("psi must be between 0 and 1.")
  ## Check that psi has proper dimensions.
  if(!identical(dim(psi),dim(x))) stop("psi has improper dimensions.")
  
  # Check log parameter.
  ## Check that log is a vector.
  if(!is.vector(log)) stop("log must be a vector.")
  ## Check that log is logical.
  if(!is.logical(log)) stop("log must be logical.")
  ## Check that log has length 1.
  if(length(log)!=1) stop("log must have length 1.")
  ## Check that log is not NA.
  if(is.na(log)) stop("log cannot be NA.")
  
  # Compute log probability density.
  d<-rowSums(stats::dbinom(x=x,size=1,prob=psi,log=TRUE))
  
  # Exponentiate if log is FALSE.
  if(!log) d<-exp(d)
  
  # Return density.
  return(d)
  
}