#' WAIC for Dirichlet-Multinomial Regression Models
#'
#' @description Computes the widely applicable information criterion (WAIC) for Dirichlet-multinomial regression models. Serves as a wrapper for [`dmreg`][dmreg()], [`dmpredict`][dmpredict()], [`ddirmult`][ddirmult()], and [`waic`][waic()] for convenient WAIC calculations. Installation of the `rstan` package is required to use this function.
#' @details For convenience, wraps the steps involved in WAIC calculations for Bayesian Dirichlet-multinomial regression models. Begins by fitting a Bayesian Dirichlet-multinomial regression model with the [`dmreg`][dmreg()] function, then generates resubstitution posterior predictions using the [`dmpredict`][dmpredict()] function. The pointwise log-likelihood is calculated with the [`ddirmult`][ddirmult()] function given the response matrix, posterior predictions, and precision parameter. WAIC is calculated from the pointwise log-likelihood using the [`waic`][waic()] function.
#' 
#' @param Y Numeric response matrix. Each record represents an observation, and each field represents a response dimension. Matrix cells contain integer counts.
#' @param X Numeric predictor matrix. Each record represents an observation, and each field represents a predictor variable. Matrix cells contain predictor values.
#' @param H Numeric vector or matrix (optional). If provided, then hierarchical effects are included in the model. Vector or matrix elements contain integer identifiers for values of hierarchical variables. If vector, then a single hierarchical variable is included, with each element representing an observation. If matrix, then each record represents an observation, and each field represents a hierarchical variable. Up to four hierarchical variables are supported (each with an arbitrary number of hierarchical levels).
#' @param ones Logical scalar. If `TRUE` (the default), then one is added to each cell of the response matrix. This avoids numerical errors which occur when distributional parameters in the model approach zero. For more information, see Harrison *et al.* (2020). If the response matrix contains no zeros, then `ones` may be set to `FALSE`.
#' @param method Numeric scalar. Options are `1` or `2`, representing the alternative WAIC bias correction formulas (*p*WAIC1 and *p*WAIC2, respectively) described in Gelman *et al.* (2014). As recommended by Gelman *et al.* (2014), the default `method` (`2`) uses the *p*WAIC2 bias correction formula.
#' @param priors Named numeric vector. Elements represent the prior values of their respective named parameters. When predictors are centered and scaled, the defaults generally represent weakly informative priors. Regression coefficients (`B`) and the precision parameter (`theta`) receive normal priors (with standard normal as the default). If hierarchical variables (argument `H`) are provided, then the common variances receive inverse-gamma priors (with default `alpha` and `beta` parameters of `0.01`).
#' @param control Named list of parameters which control the behavior of the Stan sampler. Passed to the `control` argument of the `rstan::sampling` function.
#' @param ... Additional arguments passed to the `rstan::sampling` function.
#' @returns Returns numeric scalar of the widely applicable information criterion.
#' @seealso
#' [`dmreg`][dmreg()] for fitting Dirichlet-multinomial regression models. \cr \cr
#' [`dmpredict`][dmpredict()] for generating predictions from Dirichlet-multinomial regression models. \cr \cr
#' [`ddirmult`][ddirmult()] for probability mass function of the Dirichlet-multinomial distribution. \cr \cr
#' [`waic`][waic()] for generic function to compute widely applicable information criterion.
#' @references Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M, Brubaker M, Guo J, Li P, and Riddell A. 2017. Stan: A probabilistic programming language. *Journal of Statistical Software*, 76: 1-32. DOI: 10.18637/jss.v076.i01 \cr \cr
#' Gelman A, Hwang J, and Vehtari A. 2014. Understanding predictive information criteria for Bayesian models. *Statistics and Computing*, 24(6): 997-1016. DOI: 10.1007/s11222-013-9416-2 \cr \cr
#' Goodwin KB, Hutchinson JD, and Gompert Z. 2022. Spatiotemporal and ontogenetic variation, microbial selection, and predicted *Bd*-inhibitory function in the skin-associated microbiome of a Rocky Mountain amphibian. *Frontiers in Microbiology*, 13: 1020329. DOI: 10.3389/fmicb.2022.1020329 \cr \cr
#' Harrison JG, Calder WJ, Shastry V, and Buerkle CA. Dirichlet-multinomial modelling outperforms alternatives for analysis of microbiome and other ecological count data. *Molecular Ecology Resources*, 20(2): 481-497. DOI: 10.1111/1755-0998.13128 \cr \cr
#' Watanabe S. 2010. Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11(116): 3571-3594.
#' @examplesIf interactive()
#' # Define example data file path.
#' path<-system.file("extdata",
#'                   "example_regression_data.rds",
#'                   package="LocaTT",
#'                   mustWork=TRUE)
#' 
#' # Read in example regression data.
#' data<-readRDS(file=path)
#' 
#' # Compute WAIC for Dirichlet-multinomial regression.
#' out<-dmWAIC(Y=data$Y,X=data$X,H=data$H)
#' @export
dmWAIC<-function(Y,X,H,ones=TRUE,method=2,priors=c(B.mu=0,B.sd=1,theta.mu=0,theta.sd=1,sigma2.alpha=0.01,sigma2.beta=0.01),control=list(adapt_delta=0.95,max_treedepth=20),...){
  
  # Ensure that the rstan package is installed.
  if(!requireNamespace(package="rstan",quietly=TRUE)){
    stop("Please install package 'rstan' to use this function.")
  }
  
  # Check Y input data.
  ## Check for matrix format.
  if(!is.matrix(Y)) stop("Y must be a matrix.")
  ## Check for at least two columns.
  if(ncol(Y) < 2) stop("Y must have at least two columns.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(Y))) stop("Y cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(Y)) stop("Y must be numeric integer.")
  ## Check that matrix is integer.
  if(any(Y %% 1 != 0)) stop("Y must be numeric integer.")
  ## Check that values are positive.
  if(any(Y < 0)) stop("Y must be positive.")
  
  # Check ones input data.
  ## Check that ones is a vector.
  if(!is.vector(ones)) stop("ones must be a vector.")
  ## Check that ones is logical.
  if(!is.logical(ones)) stop("ones must be logical.")
  ## Check that ones has length 1.
  if(length(ones)!=1) stop("ones must have length 1.")
  ## Check that ones is not NA.
  if(is.na(ones)) stop("ones cannot be NA.")
  
  # Add ones.
  if(ones) Y<-Y+1
  
  # If there are no hierarchical variables.
  if(missing(H)){
    
    # Fit regression.
    fit<-dmreg(Y=Y,X=X,ones=FALSE,priors=priors,control=control,...)
    
    # Generate predictions.
    P<-dmpredict(X=X,fit=fit)
    
  }else{ # If there are hierarchical variables.
    
    # Fit regression.
    fit<-dmreg(Y=Y,X=X,H=H,ones=FALSE,priors=priors,control=control,...)
    
    # Generate predictions.
    P<-dmpredict(X=X,H=H,fit=fit)
    
  }
  
  # Retrieve posterior samples.
  MCMC<-as.matrix(fit)
  
  # Retrieve precision parameter.
  theta<-MCMC[,"theta"]
  
  # Retrieve number of observations.
  N<-nrow(X)
  
  # Retrieve number of MCMC samples.
  S<-nrow(MCMC)
  
  # Define function for log probability density.
  logprob<-function(i){
    # Format response data.
    x<-Y[rep(x=i,times=S),]
    # Retrieve predictions.
    p<-P[[i]]
    # Compute log probability density.
    d<-ddirmult(x=x,p=p,theta=theta,log=TRUE)
    # Return log probability density.
    return(d)
  }
  
  # Compute log probability densities.
  loglik<-sapply(X=1:N,FUN=logprob)
  
  # Compute WAIC.
  out<-waic(loglik=loglik,method=method)
  
  # Return WAIC.
  return(out)
  
}