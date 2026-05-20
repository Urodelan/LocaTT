#' WAIC for Multivariate Logistic Regression Models
#'
#' @description Computes the widely applicable information criterion (WAIC) for multivariate logistic regression models. Serves as a wrapper for [`mlreg`][mlreg()], [`mlpredict`][mlpredict()], [`djsdm`][djsdm()], and [`waic`][waic()] for convenient WAIC calculations. Installation of the `rstan` package is required to use this function.
#' @details For convenience, wraps the steps involved in WAIC calculations for Bayesian multivariate logistic regression models. Begins by fitting a Bayesian multivariate logistic regression model with the [`mlreg`][mlreg()] function, then generates resubstitution posterior predictions using the [`mlpredict`][mlpredict()] function. The pointwise log-likelihood is calculated with the [`djsdm`][djsdm()] function given the response matrix and posterior predictions. WAIC is calculated from the pointwise log-likelihood using the [`waic`][waic()] function. Because [`djsdm`][djsdm()] does not consider residual correlations in density calculations, species interactions do not contribute to WAIC (*i.e.*, response dimensions are independent).
#' 
#' @param Y Numeric response matrix. Each record represents an observation, and each field represents a response dimension. Matrix cells contain binary values (*i.e.*, `0` or `1`).
#' @param X Numeric predictor matrix. Each record represents an observation, and each field represents a predictor variable. Matrix cells contain predictor values.
#' @param multivariate Logical scalar. If `TRUE` (the default), then fits a multivariate logistic regression. If `FALSE`, then fits stacked univariate logistic regressions.
#' @param method Numeric scalar. Options are `1` or `2`, representing the alternative WAIC bias correction formulas (*p*WAIC1 and *p*WAIC2, respectively) described in Gelman *et al.* (2014). As recommended by Gelman *et al.* (2014), the default `method` (`2`) uses the *p*WAIC2 bias correction formula.
#' @param priors Named numeric vector. Elements represent the prior values of their respective named parameters. When predictors are centered and scaled, the defaults generally represent weakly informative priors. Regression coefficients (`B`) receive normal priors (with standard normal as the default). The residual correlation matrix (`R`) receives an LKJ prior (with default shape parameter of `1`).
#' @param iter Numeric scalar. Integer value specifying the number of iterations for each chain (including warmup). The default is `20000`. Passed to the `iter` argument of the `rstan::sampling` function.
#' @param thin Numeric scalar. Integer value specifying the thinning interval. The default is `20`. Passed to the `thin` argument of the `rstan::sampling` function.
#' @param control Named list of parameters which control the behavior of the Stan sampler. Passed to the `control` argument of the `rstan::sampling` function.
#' @param ... Additional arguments passed to the `rstan::sampling` function.
#' @returns Returns numeric scalar of the widely applicable information criterion.
#' @seealso
#' [`mlreg`][mlreg()] for fitting multivariate logistic regression models. \cr \cr
#' [`mlpredict`][mlpredict()] for generating predictions from multivariate logistic regression models. \cr \cr
#' [`djsdm`][djsdm()] for probability mass function of a joint species distribution model. \cr \cr
#' [`waic`][waic()] for generic function to compute widely applicable information criterion.
#' @references Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M, Brubaker M, Guo J, Li P, and Riddell A. 2017. Stan: A probabilistic programming language. *Journal of Statistical Software*, 76: 1-32. DOI: 10.18637/jss.v076.i01 \cr \cr
#' Gelman A, Hwang J, and Vehtari A. 2014. Understanding predictive information criteria for Bayesian models. *Statistics and Computing*, 24(6): 997-1016. DOI: 10.1007/s11222-013-9416-2 \cr \cr
#' O'Brien SM, and Dunson DB. 2004. Bayesian multivariate logistic regression. *Biometrics*, 60: 739-746. DOI: 10.1111/j.0006-341X.2004.00224.x \cr \cr
#' Ovaskainen O, Hottola J, and Siitonen J. 2010. Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. *Ecology*, 91(9): 2514-2521. DOI: 10.1890/10-0173.1 \cr \cr
#' Watanabe S. 2010. Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11(116): 3571-3594.
#' @examplesIf interactive()
#' # Define example data file path.
#' path<-system.file("extdata",
#'                   "example_mvlogistic_data.rds",
#'                   package="LocaTT",
#'                   mustWork=TRUE)
#' 
#' # Read in example regression data.
#' data<-readRDS(file=path)
#' 
#' # Compute WAIC for multivariate logistic regression.
#' out<-mlWAIC(Y=data$Y,X=data$X)
#' @export
mlWAIC<-function(Y,X,multivariate=TRUE,method=2,priors=c(B.mu=0,B.sd=1,lkj=1),iter=20000,thin=20,control=list(adapt_delta=0.99,max_treedepth=20,stepsize=0.01),...){
  
  # Ensure that the rstan package is installed.
  if(!requireNamespace(package="rstan",quietly=TRUE)){
    stop("Please install package 'rstan' to use this function.")
  }
  
  # Fit regression.
  fit<-mlreg(Y=Y,X=X,multivariate=multivariate,
             priors=priors,iter=iter,thin=thin,
             control=control,...)
  
  # Generate predictions.
  P<-mlpredict(X=X,fit=fit)
  
  # Retrieve posterior samples.
  MCMC<-as.matrix(fit)
  
  # Retrieve number of observations.
  N<-nrow(X)
  
  # Retrieve number of MCMC samples.
  S<-nrow(MCMC)
  
  # Define function for log probability density.
  logprob<-function(i){
    # Format response data.
    x<-Y[rep(x=i,times=S),]
    # If single dimension.
    if(ncol(Y)==1){
      # Convert from vector to matrix.
      x<-matrix(data=x,ncol=1)
    }
    # Retrieve predictions.
    psi<-P[[i]]
    # Compute log probability density.
    d<-djsdm(x=x,psi=psi,log=TRUE)
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