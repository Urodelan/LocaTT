#' Bray-Curtis Dissimilarity (Proportion Form)
#'
#' @description Compute Bray-Curtis dissimilarity from proportional abundances.
#' @details Calculates Bray-Curtis dissimilarity from proportional abundances using the formula [`sum`][sum()]`(`[`abs`][abs()]`(p1-p2))/2`. This is equivalent to the Bray-Curtis dissimilarity formula in the `vegdist` function of the `vegan` package when both communities have the same total counts (*e.g.*, rarefied counts). This function is primarily intended to provide a method to compute Bray-Curtis dissimilarity from the posterior predictions of Dirichlet-multinomial regression models, which generate predictions of proportional abundances.
#' 
#' The dimensions of `p1` and `p2` must match. If `p1` and `p2` are matrices, then each record represents paired replicates for the communities (*e.g.*, predictions from the same posterior draw). The elements (if `p1` and `p2` are vectors) or fields (if `p1` and `p2` are matrices) represent dimensions (*e.g.*, taxa or species). If `p1` and `p2` are vectors, then a numeric scalar is returned. If `p1` and `p2` are matrices, then a numeric vector is returned whose elements correspond to the paired records of `p1` and `p2`.
#' @param p1 Numeric vector or matrix of proportional abundances for first community. See details.
#' @param p2 Numeric vector or matrix of proportional abundances for second community. See details.
#' @returns Numeric scalar or vector of Bray-Curtis dissimilarity values.
#' @seealso
#' [`dmreg`][dmreg()] for fitting Dirichlet-multinomial regression models. \cr \cr
#' [`dmpredict`][dmpredict()] for generating predictions from Dirichlet-multinomial regression models. \cr \cr
#' [`diversity`][diversity()] for computing Hill diversity from proportional abundances. \cr \cr
#' [`richness`][richness()] for computing species richness from occupancy probabilities.
#' @references Bray JR, and Curtis JT. 1957. An ordination of the upland forest communities of southern Wisconsin. *Ecological Monographs*, 27(4): 325-349. DOI: 10.2307/1942268 \cr \cr
#' Legendre P, and Legendre L. 2012. *Numerical Ecology: Third Edition*. Elsevier. \cr \cr
#' Odum EP. 1950. Bird populations of the Highlands (North Carolina) Plateau in relation to plant succession and avian invasion. *Ecology*, 31(4): 587-605. DOI: 10.2307/1931577
#' @examples
#' # Compute Bray-Curtis dissimilarity.
#' dissimilarity(p1=c(0.15,0.25,0.4,0.2),
#'               p2=c(0.25,0.35,0.1,0.3))
#' @export
dissimilarity<-function(p1,p2){
  
  # Check p1 parameter.
  ## If p1 is a vector.
  if(is.vector(p1)){
    ## Convert vector to matrix.
    p1<-matrix(p1,nrow=1,byrow=FALSE)
  }
  ## Check that p1 is a matrix.
  if(!is.matrix(p1)) stop("p1 must be a vector or matrix.")
  ## Check that p1 is numeric.
  if(!is.numeric(p1)) stop("p1 must be numeric.")
  ## Check that p1 is positive.
  if(any(p1 < 0)) stop("p1 must be positive.")
  
  # Check p2 parameter.
  ## If p2 is a vector.
  if(is.vector(p2)){
    ## Convert vector to matrix.
    p2<-matrix(p2,nrow=1,byrow=FALSE)
  }
  ## Check that p2 is a matrix.
  if(!is.matrix(p2)) stop("p2 must be a vector or matrix.")
  ## Check that p2 is numeric.
  if(!is.numeric(p2)) stop("p2 must be numeric.")
  ## Check that p2 is positive.
  if(any(p2 < 0)) stop("p2 must be positive.")
  
  # Check for matching dimensions between p1 and p2.
  if(!identical(dim(p1),dim(p2))) stop("Dimension mistmach.")
  
  # Normalize p1 to ensure simplex.
  p1<-normalize(x=p1)
  
  # Normalize p2 to ensure simplex.
  p2<-normalize(x=p2)
  
  # Define internal dissimilarity function.
  fun<-function(p1,p2) sum(abs(p1-p2))/2
  
  # Vectorize dissimilarity calculations.
  bc<-sapply(X=1:nrow(p1),FUN=function(i) fun(p1=p1[i,],p2=p2[i,]))
  
  # Return dissimilarity.
  return(bc)
  
}