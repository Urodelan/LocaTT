#' Species Richness
#'
#' @description Compute species richness from occupancy probabilities.
#' @details Calculates species richness from occupancy probabilities. Given a vector of species occupancy probabilities, computes the expected number of species as the sum of the probabilities. If given a matrix of species occupancy probabilities (where each record represents a community), computes the expected number of species as the row sums.
#' @param psi Numeric vector or matrix of occupancy probabilities. If vector, then species richness is computed for the vector of probabilities (and a scalar is returned). If matrix, then species richness is computed independently for each record (and a vector is returned).
#' @returns Numeric scalar or vector of species richness values.
#' @seealso
#' [`diversity`][diversity()] for computing Hill diversity from proportional abundances. \cr \cr
#' [`dissimilarity`][dissimilarity()] for computing Bray-Curtis dissimilarity from proportional abundances.
#' @examples
#' # Compute species richness.
#' richness(psi=c(0.506,0.825,0.135,0.683))
#' @export
richness<-function(psi){
  
  # Check psi parameter.
  ## If psi is a vector.
  if(is.vector(psi)){
    ## Convert vector to matrix.
    psi<-matrix(psi,nrow=1,byrow=FALSE)
  }
  ## Check that psi is a matrix.
  if(!is.matrix(psi)) stop("psi must be a vector or matrix.")
  ## Check that psi is numeric.
  if(!is.numeric(psi)) stop("psi must be numeric.")
  ## Check that psi is between 0 and 1.
  if(any((psi < 0) | (psi > 1),na.rm=TRUE)) stop("psi must be between 0 and 1.")
  
  # Compute richness.
  spp<-rowSums(psi)
  
  # Return richness.
  return(spp)
  
}