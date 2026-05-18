#' Predictions for Multivariate Logistic Regression Models
#'
#' @description Generate predictions for multivariate logistic regression models.
#' @details Generates posterior predictions for multivariate logistic regression models fit with the [`mlreg`][mlreg()] function. Returns a list where each element contains a matrix of posterior predictions for the respective record of `X`. Field names for the element matrices can optionally be provided with the `names` argument.
#' @param X Numeric predictor matrix. Predictions are made for each record. Each field represents a predictor variable, and the predictor variables must match (in order) those used to fit the [`mlreg`][mlreg()] model. Matrix cells contain predictor values. Element names in the returned list are taken from the row names of `X`.
#' @param fit A `stanfit` object returned from the [`mlreg`][mlreg()] function. The fitted multivariate logistic regression model.
#' @param names Vector (optional). If provided, then field names in the matrices of the returned list will receive these values. If omitted, then the matrices in the returned list will lack field names.
#' @returns A list whose elements contain numeric matrices of posterior predictions. Within the list, one element is returned for each record of `X`. Element names are taken from the row names of `X`.
#' @seealso
#' [`mlreg`][mlreg()] for fitting multivariate logistic regression models. \cr \cr
#' [`mlWAIC`][mlWAIC()] for computing widely applicable information criteria for multivariate logistic regression models.
#' @examples
#' # Define example data file path.
#' path<-system.file("extdata",
#'                   "example_mvlogistic_data.rds",
#'                   package="LocaTT",
#'                   mustWork=TRUE)
#' 
#' # Read in example regression data.
#' data<-readRDS(file=path)
#' 
#' # Predict with fitted multivariate logistic regression.
#' out<-mlpredict(X=data$X,fit=data$fit,names=colnames(data$Y))
#' @export
mlpredict<-function(X,fit,names){
  
  # Check X input data.
  ## Check for matrix format.
  if(!is.matrix(X)) stop("X must be a matrix.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(X))) stop("X cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(X)) stop("X must be numeric.")
  
  # Check that model is a stanfit object.
  if(!inherits(x=fit,what="stanfit")) stop("Model must be a 'stanfit' object.")
  
  # Check that model is a mlreg fit.
  if(fit@model_name!="mlreg") stop("Model must be a 'mlreg' fit.")
  
  # Retrieve posterior samples.
  MCMC<-as.matrix(fit)
  
  # Subset to regression coefficients.
  MCMC.B<-MCMC[,grepl(pattern="^B\\[",x=colnames(MCMC)),drop=FALSE]
  
  # Retrieve number of predictors.
  K<-ncol(X)
  
  # Retrieve number of dimensions.
  D<-ncol(MCMC.B)/K
  
  # Retrieve number of posterior samples.
  S<-nrow(MCMC.B)
  
  # Retrieve number of observations.
  N<-nrow(X)
  
  # Verify number of predictor variables.
  ## Removing leading characters.
  predictor<-sub(pattern="^B\\[[[:digit:]]+,",replacement="",x=colnames(MCMC.B))
  ## Remove trailing characters.
  predictor<-sub(pattern="\\]$",replacement="",x=predictor)
  ## Format values as numeric.
  predictor<-as.numeric(predictor)
  ## Get number of predictor variables.
  predictor<-max(predictor)
  ## Check number of predictor variables.
  if(K!=predictor){
    stop("Improper number of predictor variables.")
  }
  
  # Verify number of response variables.
  ## Removing leading characters.
  response<-sub(pattern="^B\\[",replacement="",x=colnames(MCMC.B))
  ## Remove trailing characters.
  response<-sub(pattern="\\,[[:digit:]]+]$",replacement="",x=response)
  ## Format values as numeric.
  response<-as.numeric(response)
  ## Get number of response variables.
  response<-max(response)
  ## Check number of response variables.
  if(D!=response){
    stop("Improper number of response variables.")
  }
  
  # If column names are provided.
  if(!missing(names)){
    
    # Check that column names is a vector.
    if(!is.vector(names)) stop("Names must be a vector.")
    
    # Check that column names has proper length.
    if(length(names)!=D) stop("Names has improper length.")
    
  }
  
  # Create sample-wise storage list.
  store<-vector(mode="list",length=S)
  
  # Loop through each posterior sample.
  for(s in 1:S){
    
    # Retrieve regression coefficient matrix.
    B<-matrix(data=MCMC.B[s,],nrow=D,ncol=K,byrow=FALSE)
    
    # Compute linear combination.
    eta<-X %*% t(B)
    
    # Store linear combination.
    store[[s]]<-eta
    
  }
  
  # Create observation-wise storage list.
  P<-vector(mode="list",length=N)
  
  # Loop through each observation.
  for(i in 1:N){
    
    # Retrieve posterior predictions.
    P[[i]]<-t(sapply(X=store,FUN=function(x) x[i,]))
    
    # If single dimension.
    if(D==1){
      
      # Transpose again.
      P[[i]]<-t(P[[i]])
      
      # Clear row names.
      row.names(P[[i]])<-NULL
      
    }
    
    # Apply logit link function.
    P[[i]]<-stats::plogis(P[[i]])
    
  }
  
  # If column names are provided.
  if(!missing(names)){
    
    # Define function for setting column names.
    name.cols<-function(x){
      # Set column names.
      colnames(x)<-names
      # Return named object.
      return(x)
    }
    
    # Set column names for each element.
    P<-sapply(X=P,FUN=name.cols,simplify=FALSE)
    
  }
  
  # Provide element names.
  names(P)<-row.names(X)
  
  # Return posterior predictions.
  return(P)
  
}