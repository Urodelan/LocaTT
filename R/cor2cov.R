#' Convert Correlation to Covariance Matrix
#'
#' @description Derives covariance matrix from correlation matrix and standard deviation vector.
#' @details Given correlation matrix `R` and standard deviation vector `sd`, performs the operation [`diag`][diag()]`(sd)` [`%*%`][matmult()] `R` [`%*%`][matmult()] [`diag`][diag()]`(sd)` to derive the corresponding covariance matrix. This is a counterpart to [`stats::cov2cor`][stats::cov2cor()], which scales a covariance matrix into the corresponding correlation matrix.
#' @param sd Numeric vector of standard deviations.
#' @param R Numeric correlation matrix.
#' @returns Returns a numeric covariance matrix.
#' @seealso
#' [`stats::cov2cor`][stats::cov2cor()] for scaling a covariance matrix into the corresponding correlation matrix.
#' @examples
#' # Define standard deviation vector.
#' sd<-c(9.655,1.157,1.128,2.925)
#' 
#' # Define correlation matrix.
#' R<-matrix(data=c(1.000,-0.80,0.64,-0.512,
#'                  -0.800,1.00,-0.80,0.640,
#'                  0.640,-0.80,1.00,-0.800,
#'                  -0.512,0.64,-0.80,1.000),
#'            ncol=4,byrow=TRUE)
#' 
#' # Derive covariance matrix.
#' cor2cov(sd=sd,R=R)
#' @export
cor2cov<-function(sd,R){
  
  # Check sd parameter.
  ## Check that sd is a vector.
  if(!is.vector(sd)) stop("sd must be a vector.")
  ## Check that sd is numeric.
  if(!is.numeric(sd)) stop("sd must be numeric.")
  ## Check that sd is positive.
  if(any(sd < 0,na.rm=TRUE)) stop("sd must be positive.")
  
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
  if(ncol(R)!=length(sd)) stop("R has improper dimensions.")
  
  # Compute covariance.
  Sigma<-diag(sd) %*% R %*% diag(sd)
  
  # Return covariance.
  return(Sigma)
  
}