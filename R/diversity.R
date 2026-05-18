#' Hill Diversity
#'
#' @description Compute Hill diversity from proportional abundances.
#' @details Calculates Hill diversity from proportional abundances as defined in Hill (1973), which provides a unifying theory for ecological diversity indices. When `alpha = 0`, Hill diversity is equal to species richness. When `alpha = 1`, Hill diversity is equal to the exponentiated Shannon's entropy. When `alpha = 2` (the default), Hill diversity is equal to the inverse of Simpson's index. For any value of `alpha`, the Hill diversity of a community with uniform proportional abundances is equal to species richness. Hill diversity represents the effective number of species.
#' @param p Numeric vector or matrix of proportional abundances. If vector, then Hill diversity is computed for the vector of proportions (and a scalar is returned). If matrix, then Hill diversity is computed independently for each record (and a vector is returned).
#' @param alpha Numeric scalar or vector. Continuous positive alpha parameter of the Hill diversity formula. If scalar, then `alpha` is recycled for each record of matrix `p`. If vector, then each element of `alpha` is applied to the corresponding record of matrix `p`. With the default of `alpha = 2`, Hill diversity is equal to the inverse Simpson index.
#' @returns Numeric scalar or vector of Hill diversity values.
#' @seealso
#' [`dissimilarity`][dissimilarity()] for computing Bray-Curtis dissimilarity from proportional abundances. \cr \cr
#' [`richness`][richness()] for computing species richness from occupancy probabilities.
#' @references Hill MO. 1973. Diversity and evenness: A unifying notation and its consequences. *Ecology*, 54(2): 427-432. DOI: 10.2307/1934352
#' @examples
#' # Compute Hill diversity.
#' diversity(p=c(0.15,0.25,0.4,0.2))
#' @export
diversity<-function(p,alpha=2){
  
  # Check p parameter.
  ## If p is a vector.
  if(is.vector(p)){
    ## Convert vector to matrix.
    p<-matrix(p,nrow=1,byrow=FALSE)
  }
  ## Check that p is a matrix.
  if(!is.matrix(p)) stop("p must be a vector or matrix.")
  ## Check that p is numeric.
  if(!is.numeric(p)) stop("p must be numeric.")
  ## Check that p is positive.
  if(any(p < 0)) stop("p must be positive.")
  
  # Check alpha parameter.
  ## Check that alpha is a vector.
  if(!is.vector(alpha)) stop("alpha must be a vector.")
  ## If alpha has length one.
  if(length(alpha)==1){
    ## Match length with p.
    alpha<-rep(x=alpha,times=nrow(p))
  }
  ## Check that alpha has proper length.
  if(length(alpha)!=nrow(p)) stop("alpha has improper length.")
  ## Check that alpha is numeric.
  if(!is.numeric(alpha)) stop("alpha must be numeric.")
  ## Check that alpha is positive.
  if(any(alpha < 0)) stop("alpha must be positive.")
  
  # Normalize p to ensure simplex.
  p<-normalize(x=p)
  
  # Define internal diversity function.
  fun<-function(p,alpha){
    
    # If alpha is one.
    if(alpha==1){
      
      # Compute weighted geometric mean.
      div<-exp(-sum(p*log(p)))
      
    }else{ # If alpha is not one.
      
      # Compute hill number.
      div<-sum(p^alpha)^(1/(1-alpha))
      
    }
    
    # Return diversity.
    return(div)
    
  }
  
  # Vectorize diversity calculations.
  hill<-sapply(X=1:nrow(p),FUN=function(i) fun(p=p[i,],alpha=alpha[i]))
  
  # Return diversity.
  return(hill)
  
}