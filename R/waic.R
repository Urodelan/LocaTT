#' Widely Applicable Information Criterion
#'
#' @description Generic function calculating widely applicable information criterion (WAIC) from the pointwise log-likelihood.
#' @details Given the pointwise log-likelihood, calculates WAIC (Watanabe 2010) using the formulas described in Gelman *et al.* (2014). The expected log pointwise predictive density (elppd) is estimated as the log pointwise predictive density (lppd) adjusted by a bias correction (either *p*WAIC1 or *p*WAIC2). To reflect the deviance scale, WAIC is defined as the elppd times negative two. As recommended by Gelman *et al.* (2014), *p*WAIC2 is used as the default bias correction (`method` = `2`). See Gelman *et al.* (2014) for details.
#' @param loglik Numeric matrix of the pointwise log-likelihood. Each record represents a Markov chain Monte Carlo (MCMC) sample, and each field represents an observation.
#' @param method Numeric scalar. Options are `1` or `2`, representing the alternative WAIC bias correction formulas (*p*WAIC1 and *p*WAIC2, respectively) described in Gelman *et al.* (2014). As recommended by Gelman *et al.* (2014), the default `method` (`2`) uses the *p*WAIC2 bias correction formula.
#' @returns Returns numeric scalar of the widely applicable information criterion.
#' @seealso
#' [`dmWAIC`][dmWAIC()] for computing widely applicable information criteria for Dirichlet-multinomial regression models. \cr \cr
#' [`mlWAIC`][mlWAIC()] for computing widely applicable information criteria for multivariate logistic regression models.
#' @references Gelman A, Hwang J, and Vehtari A. 2014. Understanding predictive information criteria for Bayesian models. *Statistics and Computing*, 24(6): 997-1016. DOI: 10.1007/s11222-013-9416-2 \cr \cr
#' Watanabe S. 2010. Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11(116): 3571-3594.
#' @examples
#' # Define example data file path.
#' path<-system.file("extdata",
#'                   "example_regression_data.rds",
#'                   package="LocaTT",
#'                   mustWork=TRUE)
#' 
#' # Read in example regression data.
#' data<-readRDS(file=path)
#' 
#' # Compute WAIC from pointwise log-likelihood.
#' out<-waic(loglik=data$loglik)
#' @export
waic<-function(loglik,method=2){
  
  # Check loglik input.
  ## Check that loglik is a matrix.
  if(!is.matrix(loglik)) stop("loglik must be a matrix.")
  ## Check that loglik is numeric.
  if(!is.numeric(loglik)) stop("loglik must be numeric.")
  ## Check that loglik does not contain NAs.
  if(any(is.na(loglik))) stop("loglik cannot contain NAs.")
  
  # Check method input.
  ## Check that method is a vector.
  if(!is.vector(method)) stop("method must be a vector.")
  ## Check that method is numeric.
  if(!is.numeric(method)) stop("method must be numeric.")
  ## Check that method has length one.
  if(length(method)!=1) stop("method must have length one.")
  ## Check that method is either one or two.
  if(!(method %in% 1:2)) stop("method must be 1 or 2.")
  
  # Compute log pointwise predictive density.
  lppd<-sum(log(colMeans(exp(loglik))))
  
  # If method is one.
  if(method==1){
    
    # Compute bias correction.
    pwaic<-2*sum(log(colMeans(exp(loglik)))-colMeans(loglik))
    
  }else{ # If method is two.
    
    # Compute bias correction.
    pwaic<-sum(apply(X=loglik,MARGIN=2,FUN=stats::var))
    
  }
  
  # Compute expected log pointwise predictive density.
  elppd<-lppd-pwaic
  
  # Define WAIC as elppd times negative two.
  out<-elppd*-2
  
  # Return WAIC.
  return(out)
  
}