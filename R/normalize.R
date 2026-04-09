#' Normalize a Vector or Matrix
#'
#' @description Normalizes a vector or each record of a matrix into a simplex.
#' @details Returns a vector or matrix whose elements (if vector) or records (if matrix) are computed as `x/`[`sum`][sum()]`(x)`. This normalizes a vector or matrix record into a set of proportions which sum to one (*i.e.*, a simplex). If a matrix is provided for the `x` argument, then normalization is performed independently for each record.
#' @param x Numeric vector or matrix. If vector, then the vector will be normalized to sum to one. If matrix, then each record will be normalized to sum to one (and a matrix returned).
#' @returns A numeric vector or matrix whose elements (if vector) or records (if matrix) sum to one.
#' @seealso
#' [`softmax`][softmax()] for the softmax function.
#' @examples
#' # Normalize vector.
#' normalize(x=c(3,1,5,7))
#' @export
normalize<-function(x){
  
  # Check x parameter.
  ## Check that x is a vector or matrix.
  if(!(is.vector(x) | is.matrix(x))) stop("x must be a vector or matrix.")
  ## Check that x is numeric.
  if(!is.numeric(x)) stop("x must be numeric.")
  
  # Define internal normalization function.
  fun<-function(x) x/sum(x)
  
  # If x is a vector.
  if(is.vector(x)){
    
    # Compute vector normalization.
    n<-fun(x=x)
    
  }else{ # If x is a matrix.
    
    # Compute matrix normalization.
    n<-t(apply(X=x,MARGIN=1,FUN=fun))
    
  }
  
  # Return normalization.
  return(n)
  
}