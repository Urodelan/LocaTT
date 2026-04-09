#' The Softmax
#'
#' @description Applies the softmax to a vector or each record of a matrix.
#' @details Returns a vector or matrix whose elements (if vector) or records (if matrix) are computed as [`exp`][exp()]`(x)/`[`sum`][sum()]`(`[`exp`][exp()]`(x))`. The softmax converts a vector or matrix record into a set of proportions which sum to one. If a matrix is provided for the `x` argument, then the softmax is applied independently for each record.
#' @param x Numeric vector or matrix. If vector, then the softmax of the vector will be returned. If matrix, then the softmax will be applied independently to each record (and a matrix returned).
#' @returns A numeric vector or matrix whose elements (if vector) or records (if matrix) sum to one.
#' @seealso
#' [`normalize`][normalize()] for vector or matrix normalization.
#' @examples
#' # Perform softmax on vector.
#' softmax(x=c(-0.25,0.75,1.5,0))
#' @export
softmax<-function(x){
  
  # Check x parameter.
  ## Check that x is a vector or matrix.
  if(!(is.vector(x) | is.matrix(x))) stop("x must be a vector or matrix.")
  ## Check that x is numeric.
  if(!is.numeric(x)) stop("x must be numeric.")
  
  # Define internal softmax function.
  fun<-function(x) exp(x)/sum(exp(x))
  
  # If x is a vector.
  if(is.vector(x)){
    
    # Compute vector softmax.
    s<-fun(x=x)
    
  }else{ # If x is a matrix.
    
    # Compute matrix softmax.
    s<-t(apply(X=x,MARGIN=1,FUN=fun))
    
  }
  
  # Return softmax.
  return(s)
  
}