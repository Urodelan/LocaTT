#' Fitting Multivariate Logistic Regression Models
#'
#' @description Fit a multivariate logistic regression model. Installation of the `rstan` package is required to use this function.
#' @details Fits a multivariate logistic regression model using the `rstan` interface to Stan (Carpenter *et al.* 2017). The multivariate logistic regression follows that of Ovaskainen *et al.* 2010, where the Bernoulli marginals are reparameterized as truncated continuous latent variables (Albert & Chib 1993). The latent variables `z` receive a positive constraint when `y = 1` and a negative constraint when `y = 0`, where `z` is a linear combination of predictors with correlated standard logistic errors. Equivalently, the latent variables follow a multivariate logistic distribution with scale parameters fixed at one (O'Brien & Dunson 2004), constructed in Stan as a Gaussian copula with logistic marginals (Song 2000). A `stanfit` object of the fitted model is returned, which can be used with standard `rstan` functions to evaluate model convergence (*e.g.*, posterior trace plots, R-hat convergence diagnostics, and effective sample sizes). By default, weakly informative priors are used on the regression coefficients (`B`) and residual correlation matrix (`R`).
#' @param Y Numeric response matrix. Each record represents an observation, and each field represents a response dimension. Matrix cells contain binary values (*i.e.*, `0` or `1`).
#' @param X Numeric predictor matrix. Each record represents an observation, and each field represents a predictor variable. Matrix cells contain predictor values.
#' @param priors Named numeric vector. Elements represent the prior values of their respective named parameters. When predictors are centered and scaled, the defaults generally represent weakly informative priors. Regression coefficients (`B`) receive normal priors (with standard normal as the default). The residual correlation matrix (`R`) receives an LKJ prior (with default shape parameter of `1`).
#' @param iter Numeric scalar. Integer value specifying the number of iterations for each chain (including warmup). The default is `20000`. Passed to the `iter` argument of the `rstan::sampling` function.
#' @param thin Numeric scalar. Integer value specifying the thinning interval. The default is `20`. Passed to the `thin` argument of the `rstan::sampling` function.
#' @param control Named list of parameters which control the behavior of the Stan sampler. Passed to the `control` argument of the `rstan::sampling` function.
#' @param ... Additional arguments passed to the `rstan::sampling` function.
#' @returns Returns a `stanfit` object of the fitted multivariate logistic regression model.
#' @seealso
#' [`mlpredict`][mlpredict()] for generating predictions from multivariate logistic regression models. \cr \cr
#' [`mlWAIC`][mlWAIC()] for computing widely applicable information criteria for multivariate logistic regression models.
#' @references Albert JH, and Chib S. 1993. Bayesian analysis of binary and polychotomous response data. *Journal of the American Statistical Association*, 88(422): 669-679. \cr \cr
#' Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M, Brubaker M, Guo J, Li P, and Riddell A. 2017. Stan: A probabilistic programming language. *Journal of Statistical Software*, 76: 1-32. DOI: 10.18637/jss.v076.i01 \cr \cr
#' O'Brien SM, and Dunson DB. 2004. Bayesian multivariate logistic regression. *Biometrics*, 60: 739-746. DOI: 10.1111/j.0006-341X.2004.00224.x \cr \cr
#' Ovaskainen O, Hottola J, and Siitonen J. 2010. Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. *Ecology*, 91(9): 2514-2521. DOI: 10.1890/10-0173.1 \cr \cr
#' Song P. 2000. Multivariate dispersion models generated from Gaussian copula. *Scandinavian Journal of Statistics*, 27(2): 305-320. DOI: 10.1111/1467-9469.00191
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
#' # Fit multivariate logistic regression.
#' out<-mlreg(Y=data$Y,X=data$X)
#' @export
mlreg<-function(Y,X,priors=c(B.mu=0,B.sd=1,lkj=1),iter=20000,thin=20,control=list(adapt_delta=0.99,max_treedepth=20,stepsize=0.01),...){
  
  # Ensure that the rstan package is installed.
  if(!requireNamespace(package="rstan",quietly=TRUE)){
    stop("Please install package 'rstan' to use this function.")
  }
  
  # Check Y input data.
  ## Check for matrix format.
  if(!is.matrix(Y)) stop("Y must be a matrix.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(Y))) stop("Y cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(Y)) stop("Y must be numeric.")
  ## Check that matrix is binary.
  if(!all(Y %in% 0:1)) stop("Y must be binary.")
  
  # Check X input data.
  ## Check for matrix format.
  if(!is.matrix(X)) stop("X must be a matrix.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(X))) stop("X cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(X)) stop("X must be numeric.")
  
  # Check for matching dimensions between X and Y.
  if(nrow(Y)!=nrow(X)) stop("Dimension mismatch between X and Y.")
  
  # Check format of priors.
  ## Check that priors are provided as a vector.
  if(!is.vector(priors)) stop("Priors must be provided as named vector.")
  ## Check that priors are not a list.
  if(is.list(priors)) stop("Priors cannot be a list.")
  ## Check that priors do not cotain NAs.
  if(any(is.na(priors))) stop("Priors cannot contain NAs.")
  ## Check that priors are numeric.
  if(!is.numeric(priors)) stop("Priors must be numeric.")
  
  # Define required priors.
  check<-c("B.mu","B.sd","lkj")
  
  # If any priors are absent.
  if(!all(check %in% names(priors))){
    
    # Retrieve absent prior names.
    absent<-check[!(check %in% names(priors))]
    
    # Throw error providing absent prior names.
    stop(paste("Missing priors for variable(s):",paste(absent,collapse=", ")))
    
  }
  
  # Define model parameters.
  data<-list(
    ### Define variables.
    "N"=nrow(Y), # Number of observations.
    "D"=ncol(Y), # Number of response variables.
    "K"=ncol(X), # Number of predictors.
    "Y"=Y, # Response matrix.
    "X"=X, # Predictor matrix.
    ### Define priors.
    "B_mu"=unname(priors["B.mu"]), # Regression coefficients prior mean.
    "B_sd"=unname(priors["B.sd"]), # Regression coefficients prior standard deviation.
    "lkj"=unname(priors["lkj"]) # LKJ prior correlation.
  )
  
  # Generate initial values for regression coefficients.
  ## Initialize matrix of zeros.
  B.init<-matrix(data=0,nrow=ncol(Y),ncol=ncol(X))
  ## Compute response proportions.
  p<-colMeans(Y)
  ## Check for all ones or zeros.
  p<-ifelse(p %% 1 == 0,0.5,p)
  ## Initialize intercept at inverted values.
  B.init[,1]<-stats::qlogis(p=p)
  
  # Generate initial values for correlation matrix.
  L.init<-diag(ncol(Y))
  
  # Collect initial values in a function.
  init<-function(...) list(B=B.init,L=L.init)
  
  # Define parameters.
  pars<-c("B","R")
  
  # Define posterior check function with threshold argument
  # for Bayesian fraction of missing information.
  check_posterior<-function(fit,threshold=0.2){
    
    # Divergent transitions.
    ## Get number of divergent transitions.
    n_d<-rstan::get_num_divergent(fit)
    ## If there are divergent transitions.
    if(n_d > 0){
      ## Produce a warning.
      warning(paste0("There were ",n_d," divergent transitions after warmup. See\n",
                     "https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup\n",
                     "to find out why this is a problem and how to eliminate them."),
              call.=FALSE)
    }
    
    # Maximum tree depth.
    ## Retrieve control parameters.
    control<-fit@stan_args[[1]]$control
    ## Retrieve maximum tree depth (if present), or set to default (if absent).
    max_td<-ifelse(is.null(control),10,control$max_treedepth)
    ## Get number of transitions that exceed maximum tree depth.
    n_m<-rstan::get_num_max_treedepth(fit)
    ## If there are transitions that exceed maximum tree depth.
    if(n_m > 0){
      ## Produce a warning.
      warning(paste0("There were ",n_m,
                     " transitions after warmup that exceeded the maximum treedepth.",
                     " Increase max_treedepth above ",max_td,". See\n",
                     "https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded"),
              call.=FALSE)
    }
    
    # Bayesian fraction of missing information.
    ## Get number of chains with low Bayesian fraction of missing information.
    n_e<-sum(rstan::get_bfmi(fit) < threshold)
    ## If there are chains with low Bayesian fraction of missing information.
    if(n_e > 0){
      ## Produce a warning.
      warning(paste0("There were ",n_e,
                     " chains where the estimated Bayesian Fraction of Missing Information",
                     " was low. See\n","https://mc-stan.org/misc/warnings.html#bfmi-low"),
              call.=FALSE)
    }
    
    # If there are diagnostic issues.
    if(n_d > 0 || n_m > 0 || n_e > 0){
      # Encourage examination of pairs plot.
      warning("Examine the pairs() plot to diagnose sampling problems\n",
              call.=FALSE,noBreaks.=TRUE)
    }
    
    # Retrieve model fit as an array.
    sims<-as.array(fit)
    
    # R-hat convergence diagnostic.
    ## Compute R-hat convergence diagnostics.
    rhat<-apply(X=sims,MARGIN=3,FUN=rstan::Rhat)
    ## If any R-hat convergence diagnostics are too high.
    if(any(rhat > 1.05,na.rm=TRUE)){
      ## Produce a warning.
      warning(paste0("The largest R-hat is ",round(max(rhat,na.rm=TRUE),digits=2),
                     ", indicating chains have not mixed.\n",
                     "Running the chains for more iterations may help. See\n",
                     "https://mc-stan.org/misc/warnings.html#r-hat"),
              call.=FALSE)
    }
    
    # Bulk effective sample size.
    ## Compute bulk effective sample sizes.
    bulk_ess<-apply(X=sims,MARGIN=3,FUN=rstan::ess_bulk)
    ## If any bulk effective sample sizes are too low.
    if(any(bulk_ess < (100*ncol(sims)),na.rm=TRUE)){
      ## Produce a warning.
      warning(paste0("Bulk Effective Samples Size (ESS) is too low, ",
                     "indicating posterior means and medians may be unreliable.\n",
                     "Running the chains for more iterations may help. See\n",
                     "https://mc-stan.org/misc/warnings.html#bulk-ess"),call.=FALSE)
    }
    
    # Tail effective sample size.
    ## Compute tail effective sample sizes.
    tail_ess<-apply(X=sims,MARGIN=3,FUN=rstan::ess_tail)
    ## If any tail effective samples sizes are too low.
    if(any(tail_ess < (100*ncol(sims)),na.rm=TRUE)){
      ## Produce a warning.
      warning(paste0("Tail Effective Samples Size (ESS) is too low, indicating ",
                     "posterior variances and tail quantiles may be unreliable.\n",
                     "Running the chains for more iterations may help. See\n",
                     "https://mc-stan.org/misc/warnings.html#tail-ess"),call.=FALSE)
    }
    
  }
  
  # Define model file.
  file<-"mlreg.stan"
  
  # Define model path.
  path<-system.file("intdata",file,package="LocaTT",mustWork=TRUE)
  
  # Compile the model.
  model<-rstan::stan_model(file=path,model_name="mlreg")
  
  # Fit the model.
  suppressWarnings(
    fit<-rstan::sampling(object=model,
                         data=data,
                         pars=pars,
                         iter=iter,
                         thin=thin,
                         init=init,
                         control=control,
                         ...)
  )
  
  # Check posterior.
  check_posterior(fit=fit)
  
  # Return fit.
  return(fit)
  
}