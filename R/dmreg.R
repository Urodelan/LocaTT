#' Fitting Dirichlet-Multinomial Regression Models
#'
#' @description Fit a Bayesian Dirichlet-multinomial regression model. Both fixed and hierarchical effects are supported. Installation of the `rstan` package is required to use this function.
#' @details Fits the Bayesian Dirichlet-multinomial regression model of Goodwin *et al.* (2022) using the `rstan` interface to Stan (Carpenter *et al.* 2017). A `stanfit` object of the fitted model is returned, which can be used with standard `rstan` functions to evaluate model convergence (*e.g.*, posterior trace plots, R-hat convergence diagnostics, and effective sample sizes). The model formulation is identical to that of Goodwin *et al.* (2022), except that the hard sum-to-zero constraint on hierarchical effects was removed to preserve the prior marginal variance of the final element. Up to four hierarchical variables are supported.
#' 
#' For each observation, counts are distributed according to the Dirichlet-multinomial distribution with alpha parameters defined as the product of an expected proportions vector and an exponentiated precision parameter. The precision parameter controls the degree of overdispersion relative to the multinomial distribution. The softmax function normalizes linear predictor combinations into expected proportions. For the model to be identifiable, the regression coefficients of the final dimension are set to zero. By default, weakly informative priors are used on the regression coefficients (`B`), precision parameter (`theta`), and hierarchical variances (`sigma2`). See the supplement of Goodwin *et al.* (2022) for details.
#' @param Y Numeric response matrix. Each record represents an observation, and each field represents a response dimension. Matrix cells contain integer counts.
#' @param X Numeric predictor matrix. Each record represents an observation, and each field represents a predictor variable. Matrix cells contain predictor values.
#' @param H Numeric vector or matrix (optional). If provided, then hierarchical effects are included in the model. Vector or matrix elements contain integer identifiers for values of hierarchical variables. If vector, then a single hierarchical variable is included, with each element representing an observation. If matrix, then each record represents an observation, and each field represents a hierarchical variable. Up to four hierarchical variables are supported (each with an arbitrary number of hierarchical levels).
#' @param ones Logical scalar. If `TRUE` (the default), then one is added to each cell of the response matrix. This avoids numerical errors which occur when distributional parameters in the model approach zero. For more information, see Harrison *et al.* (2020). If the response matrix contains no zeros, then `ones` may be set to `FALSE`.
#' @param priors Named numeric vector. Elements represent the prior values of their respective named parameters. When predictors are centered and scaled, the defaults generally represent weakly informative priors. Regression coefficients (`B`) and the precision parameter (`theta`) receive normal priors (with standard normal as the default). If hierarchical variables (argument `H`) are provided, then the common variances receive inverse-gamma priors (with default `alpha` and `beta` parameters of 0.01).
#' @param control Named list of parameters which control the behavior of the Stan sampler. Passed to the `control` argument of the `rstan::sampling` function.
#' @param ... Additional arguments passed to the `rstan::sampling` function.
#' @returns Returns a `stanfit` object of the fitted Bayesian Dirichlet-multinomial regression model.
#' @seealso
#' [`dmpredict`][dmpredict()] for generating predictions from Dirichlet-multinomial regression models. \cr \cr
#' [`dmWAIC`][dmWAIC()] for computing widely applicable information criteria for Dirichlet-multinomial regression models.
#' @references Carpenter B, Gelman A, Hoffman MD, Lee D, Goodrich B, Betancourt M, Brubaker M, Guo J, Li P, and Riddell A. 2017. Stan: A probabilistic programming language. *Journal of Statistical Software*, 76: 1-32. DOI: 10.18637/jss.v076.i01 \cr \cr
#' Goodwin KB, Hutchinson JD, and Gompert Z. 2022. Spatiotemporal and ontogenetic variation, microbial selection, and predicted *Bd*-inhibitory function in the skin-associated microbiome of a Rocky Mountain amphibian. *Frontiers in Microbiology*, 13: 1020329. DOI: 10.3389/fmicb.2022.1020329 \cr \cr
#' Harrison JG, Calder WJ, Shastry V, and Buerkle CA. Dirichlet-multinomial modelling outperforms alternatives for analysis of microbiome and other ecological count data. *Molecular Ecology Resources*, 20(2): 481-497. DOI: 10.1111/1755-0998.13128
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
#' # Fit Dirichlet-multinomial regression.
#' out<-dmreg(Y=data$Y,X=data$X,H=data$H)
#' @export
dmreg<-function(Y,X,H,ones=TRUE,priors=c(B.mu=0,B.sd=1,theta.mu=0,theta.sd=1,sigma2.alpha=0.01,sigma2.beta=0.01),control=list(adapt_delta=0.95,max_treedepth=20),...){
  
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
  
  # Check X input data.
  ## Check for matrix format.
  if(!is.matrix(X)) stop("X must be a matrix.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(X))) stop("X cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(X)) stop("X must be numeric.")
  
  # Check for matching dimensions between X and Y.
  if(nrow(Y)!=nrow(X)) stop("Dimension mismatch between X and Y.")
  
  # Initialize number of hierarchical groups.
  J<-0
  
  # If hierarchical identifiers are provided.
  if(!missing(H)){
    
    # Convert vector to matrix.
    if(is.vector(H)) H<-matrix(data=H,ncol=1)
    
    # Check for matrix format.
    if(!is.matrix(H)) stop("H must be a vector or matrix.")
    
    # Check that matrix does not contain NAs.
    if(any(is.na(H))) stop("H cannot contain NAs.")
    
    # Check that matrix is numeric.
    if(!is.numeric(H)) stop("H must be numeric integer.")
    
    # Check that matrix is integer.
    if(any(H %% 1 != 0)) stop("H must be numeric integer.")
    
    # Check that values are positive.
    if(any(H <= 0)) stop("H must be positive.")
    
    # Update number of hierarchical groups.
    J<-ncol(H)
    
    # Retrieve maximum identifiers.
    maximum<-apply(X=H,MARGIN=2,FUN=max)
    
    # Check that all potential levels are included in each field.
    check<-sapply(X=1:J,FUN=function(j) all(1:maximum[j] %in% H[,j]))
    
    # If any potential hierarchical levels are missing.
    if(any(!check)){
      # Throw an error stating the problem variables.
      stop(paste("Missing levels in hierarchical variable(s):",
                 paste0(paste(which(!check),collapse=", "))))
    }
    
    # Check for matching dimensions between H and Y.
    if(nrow(Y)!=nrow(H)) stop("Dimension mismatch between H and Y.")
    
  }
  
  # Check that ones is a vector.
  if(!is.vector(ones)) stop("ones must be a vector.")
  
  # Check that ones is logical.
  if(!is.logical(ones)) stop("ones must be logical.")
  
  # Check that ones has length 1.
  if(length(ones)!=1) stop("ones must have length 1.")
  
  # Check that ones is not NA.
  if(is.na(ones)) stop("ones cannot be NA.")
  
  # Add ones.
  if(ones) Y<-Y+1
  
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
  ## Without hierarhical effects.
  fixed<-c("B.mu","B.sd","theta.mu","theta.sd")
  ## Inlcuding hierarchical effects.
  random<-c(fixed,c("sigma2.alpha","sigma2.beta"))
  
  # If no hierarchical effects.
  if(J==0){
    
    # If any priors are absent.
    if(!all(fixed %in% names(priors))){
      
      # Retrieve absent prior names.
      absent<-fixed[!(fixed %in% names(priors))]
      
      # Throw error providing absent prior names.
      stop(paste("Missing priors for variable(s):",paste(absent,collapse=", ")))
      
    }
    
    # Retrieve relevant priors.
    priors<-priors[fixed]
    
  }else{ # If including hierarchical effects.
    
    # If any priors are absent.
    if(!all(random %in% names(priors))){
      
      # Retrieve absent prior names.
      absent<-random[!(random %in% names(priors))]
      
      # Throw error providing absent prior names.
      stop(paste("Missing priors for variable(s):",paste(absent,collapse=", ")))
      
    }
    
    # Retrieve relevant priors.
    priors<-priors[random]
    
  }
  
  # Define model parameters.
  if(J==0){ # Hierarchical = 0.
    
    # Define model data.
    data<-list(
      ### Define variables.
      "N"=nrow(Y), # Number of observations.
      "D"=ncol(Y), # Number of response variables.
      "K"=ncol(X), # Number of predictors.
      "Y"=Y, # Response matrix.
      "X"=X, # Predictor matrix.
      ### Define priors.
      "B_mu"=unname(priors["B.mu"]), # Fixed effects prior mean.
      "B_sd"=unname(priors["B.sd"]), # Fixed effects prior standard deviation.
      "theta_mu"=unname(priors["theta.mu"]), # Precision parameter prior mean.
      "theta_sd"=unname(priors["theta.sd"]) # Precision parameter prior standard deviation.
    )
    
  }else if(J==1){ # Hierarchical = 1.
    
    # Retrieve hierarchical identifiers.
    ## First hierarchical variable.
    h1<-unname(H[,1])
    
    # Define model data.
    data<-list(
      ### Define variables.
      "N"=nrow(Y), # Number of observations.
      "D"=ncol(Y), # Number of response variables.
      "K"=ncol(X), # Number of predictors.
      "Y"=Y, # Response matrix.
      "X"=X, # Predictor matrix.
      "L1"=max(h1), # Number of hierarchical effects (#1).
      "h1"=h1, # Vector of hierarchical identifiers (#1).
      ### Define priors.
      "B_mu"=unname(priors["B.mu"]), # Fixed effects prior mean.
      "B_sd"=unname(priors["B.sd"]), # Fixed effects prior standard deviation.
      "theta_mu"=unname(priors["theta.mu"]), # Precision parameter prior mean.
      "theta_sd"=unname(priors["theta.sd"]), # Precision parameter prior standard deviation.
      "sigma2_alpha"=unname(priors["sigma2.alpha"]), # Hierarchical standard deviation alpha prior.
      "sigma2_beta"=unname(priors["sigma2.beta"]) # Hierarchical standard deviation beta prior.
    )
    
  }else if(J==2){ # Hierarchical = 2.
    
    # Retrieve hierarchical identifiers.
    ## First hierarchical variable.
    h1<-unname(H[,1])
    ## Second hierarchical variable.
    h2<-unname(H[,2])
    
    # Define model data.
    data<-list(
      ### Define variables.
      "N"=nrow(Y), # Number of observations.
      "D"=ncol(Y), # Number of response variables.
      "K"=ncol(X), # Number of predictors.
      "Y"=Y, # Response matrix.
      "X"=X, # Predictor matrix.
      "L1"=max(h1), # Number of hierarchical effects (#1).
      "h1"=h1, # Vector of hierarchical identifiers (#1).
      "L2"=max(h2), # Number of hierarchical effects (#2).
      "h2"=h2, # Vector of hierarchical identifiers (#2).
      ### Define priors.
      "B_mu"=unname(priors["B.mu"]), # Fixed effects prior mean.
      "B_sd"=unname(priors["B.sd"]), # Fixed effects prior standard deviation.
      "theta_mu"=unname(priors["theta.mu"]), # Precision parameter prior mean.
      "theta_sd"=unname(priors["theta.sd"]), # Precision parameter prior standard deviation.
      "sigma2_alpha"=unname(priors["sigma2.alpha"]), # Hierarchical standard deviation alpha prior.
      "sigma2_beta"=unname(priors["sigma2.beta"]) # Hierarchical standard deviation beta prior.
    )
    
  }else if(J==3){ # Hierarchical = 3.
    
    # Retrieve hierarchical identifiers.
    ## First hierarchical variable.
    h1<-unname(H[,1])
    ## Second hierarchical variable.
    h2<-unname(H[,2])
    ## Third hierarchical variable.
    h3<-unname(H[,3])
    
    # Define model data.
    data<-list(
      ### Define variables.
      "N"=nrow(Y), # Number of observations.
      "D"=ncol(Y), # Number of response variables.
      "K"=ncol(X), # Number of predictors.
      "Y"=Y, # Response matrix.
      "X"=X, # Predictor matrix.
      "L1"=max(h1), # Number of hierarchical effects (#1).
      "h1"=h1, # Vector of hierarchical identifiers (#1).
      "L2"=max(h2), # Number of hierarchical effects (#2).
      "h2"=h2, # Vector of hierarchical identifiers (#2).
      "L3"=max(h3), # Number of hierarchical effects (#3).
      "h3"=h3, # Vector of hierarchical identifiers (#3).
      ### Define priors.
      "B_mu"=unname(priors["B.mu"]), # Fixed effects prior mean.
      "B_sd"=unname(priors["B.sd"]), # Fixed effects prior standard deviation.
      "theta_mu"=unname(priors["theta.mu"]), # Precision parameter prior mean.
      "theta_sd"=unname(priors["theta.sd"]), # Precision parameter prior standard deviation.
      "sigma2_alpha"=unname(priors["sigma2.alpha"]), # Hierarchical standard deviation alpha prior.
      "sigma2_beta"=unname(priors["sigma2.beta"]) # Hierarchical standard deviation beta prior.
    )
    
  }else if(J==4){ # Hierarchical = 4.
    
    # Retrieve hierarchical identifiers.
    ## First hierarchical variable.
    h1<-unname(H[,1])
    ## Second hierarchical variable.
    h2<-unname(H[,2])
    ## Third hierarchical variable.
    h3<-unname(H[,3])
    ## Fourth hierarchical variable.
    h4<-unname(H[,4])
    
    # Define model data.
    data<-list(
      ### Define variables.
      "N"=nrow(Y), # Number of observations.
      "D"=ncol(Y), # Number of response variables.
      "K"=ncol(X), # Number of predictors.
      "Y"=Y, # Response matrix.
      "X"=X, # Predictor matrix.
      "L1"=max(h1), # Number of hierarchical effects (#1).
      "h1"=h1, # Vector of hierarchical identifiers (#1).
      "L2"=max(h2), # Number of hierarchical effects (#2).
      "h2"=h2, # Vector of hierarchical identifiers (#2).
      "L3"=max(h3), # Number of hierarchical effects (#3).
      "h3"=h3, # Vector of hierarchical identifiers (#3).
      "L4"=max(h4), # Number of hierarchical effects (#4).
      "h4"=h4, # Vector of hierarchical identifiers (#4).
      ### Define priors.
      "B_mu"=unname(priors["B.mu"]), # Fixed effects prior mean.
      "B_sd"=unname(priors["B.sd"]), # Fixed effects prior standard deviation.
      "theta_mu"=unname(priors["theta.mu"]), # Precision parameter prior mean.
      "theta_sd"=unname(priors["theta.sd"]), # Precision parameter prior standard deviation.
      "sigma2_alpha"=unname(priors["sigma2.alpha"]), # Hierarchical standard deviation alpha prior.
      "sigma2_beta"=unname(priors["sigma2.beta"]) # Hierarchical standard deviation beta prior.
    )
    
  }else{ # Hierarchical > 4.
    
    # Throw error if more than four hierarchical variables.
    stop("A maximum of 4 hierarchical variables are supported.")
    
  }
  
  # Generate initial values for fixed effects.
  ## Initialize matrix of zeros.
  B.init<-matrix(data=0,nrow=ncol(Y),ncol=ncol(X))
  ## Compute response proportions.
  p<-normalize(x=Y)
  ## Define inverse softmax function.
  inv.softmax<-function(p) log(p)-log(p[length(p)])
  ## Convert proportions to inverted values.
  inverted<-t(apply(X=p,MARGIN=1,FUN=inv.softmax))
  ## Initialize intercept at mean of inverted values.
  B.init[,1]<-colMeans(x=inverted)
  
  # Initialize precision parameter at zero.
  theta.init<-0
  
  # If hierarchical = 1+
  if(J > 0){
    
    # Generate initial values for hierarchical parameters (#1).
    ## Initialize hierarchical effects at zero.
    H1.init<-matrix(data=0,nrow=ncol(Y),ncol=max(h1))
    ## Initialize hierarchical hyperperiors at one.
    sigma2_h1.init<-rep(x=1,times=ncol(Y)-1)
    
  }
  
  # If hierarchical = 2+
  if(J > 1){
    
    # Generate initial values for hierarchical parameters (#2).
    ## Initialize hierarchical effects at zero.
    H2.init<-matrix(data=0,nrow=ncol(Y),ncol=max(h2))
    ## Initialize hierarchical hyperperiors at one.
    sigma2_h2.init<-rep(x=1,times=ncol(Y)-1)
    
  }
  
  # If hierarchical = 3+
  if(J > 2){
    
    # Generate initial values for hierarchical parameters (#3).
    ## Initialize hierarchical effects at zero.
    H3.init<-matrix(data=0,nrow=ncol(Y),ncol=max(h3))
    ## Initialize hierarchical hyperperiors at one.
    sigma2_h3.init<-rep(x=1,times=ncol(Y)-1)
    
  }
  
  # If hierarchical = 4
  if(J > 3){
    
    # Generate initial values for hierarchical parameters (#4).
    ## Initialize hierarchical effects at zero.
    H4.init<-matrix(data=0,nrow=ncol(Y),ncol=max(h4))
    ## Initialize hierarchical hyperperiors at one.
    sigma2_h4.init<-rep(x=1,times=ncol(Y)-1)
    
  }
  
  # If hierarchical = 0.
  if(J==0){
    
    # Collect initial values in a function.
    init<-function(...) list(B=B.init,theta=theta.init)
    
    # Define parameters.
    pars<-c("B","theta")
    
  }else if(J==1){ # If hierarchical = 1.
    
    # Collect initial values in a function.
    init<-function(...) list(B=B.init,theta=theta.init,
                             H1=H1.init,sigma2_h1=sigma2_h1.init)
    
    # Define parameters.
    pars<-c("B","theta",
            "H1","sigma1")
    
  }else if(J==2){ # If hierarchical = 2.
    
    # Collect initial values in a function.
    init<-function(...) list(B=B.init,theta=theta.init,
                             H1=H1.init,sigma2_h1=sigma2_h1.init,
                             H2=H2.init,sigma2_h2=sigma2_h2.init)
    
    # Define parameters.
    pars<-c("B","theta",
            "H1","sigma1",
            "H2","sigma2")
    
  }else if(J==3){ # If hierarchical = 3.
    
    # Collect initial values in a function.
    init<-function(...) list(B=B.init,theta=theta.init,
                             H1=H1.init,sigma2_h1=sigma2_h1.init,
                             H2=H2.init,sigma2_h2=sigma2_h2.init,
                             H3=H3.init,sigma2_h3=sigma2_h3.init)
    
    # Define parameters.
    pars<-c("B","theta",
            "H1","sigma1",
            "H2","sigma2",
            "H3","sigma3")
    
  }else{ # If hierarchical = 4.
    
    # Collect initial values in a function.
    init<-function(...) list(B=B.init,theta=theta.init,
                             H1=H1.init,sigma2_h1=sigma2_h1.init,
                             H2=H2.init,sigma2_h2=sigma2_h2.init,
                             H3=H3.init,sigma2_h3=sigma2_h3.init,
                             H4=H4.init,sigma2_h4=sigma2_h4.init)
    
    # Define parameters.
    pars<-c("B","theta",
            "H1","sigma1",
            "H2","sigma2",
            "H3","sigma3",
            "H4","sigma4")
    
  }
  
  # Define model file.
  file<-paste0("dmreg_h",J,".stan")
  
  # Define model path.
  path<-system.file("intdata",file,package="LocaTT",mustWork=TRUE)
  
  # Compile the model.
  model<-rstan::stan_model(file=path,model_name="dmreg")
  
  # Fit the model.
  fit<-rstan::sampling(object=model,
                       data=data,
                       pars=pars,
                       init=init,
                       control=control,
                       ...)
  
  # Return fit.
  return(fit)
  
}