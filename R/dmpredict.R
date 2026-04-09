#' Predictions for Dirichlet-Multinomial Regression Models
#'
#' @description Generate predictions for Dirichlet-multinomial regression models.
#' @details Generates posterior predictions for Dirichlet-multinomial regression models fit with the [`dmreg`][dmreg()] function. Predictions can either include or omit hierarchical effects, depending on whether argument `H` is provided. Returns a list where each element contains a matrix of posterior predictions for the respective record of `X`. Field names for the element matrices can optionally be provided with the `names` argument.
#' @param X Numeric predictor matrix. Predictions are made for each record. Each field represents a predictor variable, and the predictor variables must match (in order) those used to fit the [`dmreg`][dmreg()] model. Matrix cells contain predictor values. Element names in the returned list are taken from the row names of `X`.
#' @param H Numeric vector or matrix (optional). If provided, then hierarchical effects are included in the predictions. Vector or matrix elements contain integer identifiers for values of hierarchical variables. If vector, then a single hierarchical variable is included, with each element corresponding to a record in `X`. If matrix, then each record in `H` corresponds to a record in `X`. Each field in `H` represents a hierarchical variable, and the hierarchical variables must match (in order) those used to fit the [`dmreg`][dmreg()] model. If `H` is omitted, then hierarchical effects are omitted from the predictions (but may still have been used to fit the model).
#' @param fit A `stanfit` object returned from the [`dmreg`][dmreg()] function. The fitted Dirichlet-multinomial regression model.
#' @param names Vector (optional). If provided, then field names in the matrices of the returned list will receive these values. If omitted, then the matrices in the returned list will lack field names.
#' @returns A list whose elements contain numeric matrices of posterior predictions. Within the list, one element is returned for each record of `X`. Element names are taken from the row names of `X`.
#' @seealso
#' [`dmreg`][dmreg()] for fitting Dirichlet-multinomial regression models. \cr \cr
#' [`dmWAIC`][dmWAIC()] for computing widely applicable information criteria for Dirichlet-multinomial regression models.
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
#' # Predict with fitted Dirichlet-multinomial regression.
#' out<-dmpredict(X=data$X,fit=data$fit,names=colnames(data$Y))
#' @export
dmpredict<-function(X,H,fit,names){
  
  # Check X input data.
  ## Check for matrix format.
  if(!is.matrix(X)) stop("X must be a matrix.")
  ## Check that matrix does not contain NAs.
  if(any(is.na(X))) stop("X cannot contain NAs.")
  ## Check that matrix is numeric.
  if(!is.numeric(X)) stop("X must be numeric.")
  
  # Initialize number of hierarchical variables at zero.
  J<-0
  
  # Check H input data.
  ## If hierarchical variables are provided.
  if(!missing(H)){
    ## Convert vector to matrix.
    if(is.vector(H)) H<-matrix(data=H,ncol=1)
    ## Check for matrix format.
    if(!is.matrix(H)) stop("H must be a vector or matrix.")
    ## Check that matrix does not contain NAs.
    if(any(is.na(H))) stop("H cannot contain NAs.")
    ## Check that matrix is numeric.
    if(!is.numeric(H)) stop("H must be numeric integer.")
    ## Check that matrix is integer.
    if(any(H %% 1 != 0)) stop("H must be numeric integer.")
    ## Check that values are positive.
    if(any(H <= 0)) stop("H must be positive.")
    ## Check that dimensions match between H and X.
    if(nrow(H)!=nrow(X)) stop("Dimension mismatch between H and X.")
    ## Retrieve number of hierarchical variables.
    J<-ncol(H)
  }
  
  # Check that model is a stanfit object.
  if(!inherits(x=fit,what="stanfit")) stop("Model must be a 'stanfit' object.")
  
  # Check that model is a dmreg fit.
  if(fit@model_name!="dmreg") stop("Model must be a 'dmreg' fit.")
  
  # Retrieve posterior samples.
  MCMC<-as.matrix(fit)
  
  # Subset to regression coefficients.
  MCMC.B<-MCMC[,grepl(pattern="^B\\[",x=colnames(MCMC))]
  
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
  fixed<-sub(pattern="^B\\[[[:digit:]]+,",replacement="",x=colnames(MCMC.B))
  ## Remove trailing characters.
  fixed<-sub(pattern="\\]$",replacement="",x=fixed)
  ## Format values as numeric.
  fixed<-as.numeric(fixed)
  ## Get number of predictor variables.
  fixed<-max(fixed)
  ## Check number of predictor variables.
  if(K!=fixed){
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
  
  # If hierarchical variables are provided.
  if(J > 0){
    
    # Check if the model contains hierarchical effects.
    check<-grepl(pattern="^H[[:digit:]]+\\[",x=colnames(MCMC))
    
    # Verify that model contains hierarchical effects.
    if(!any(check)) stop("No hierarchical effects found.")
    
    # Subset to hierarchical effects.
    MCMC.H<-MCMC[,check]
    
    # Verify number of hierarchical variables.
    ## Removing leading characters.
    hierarchical<-sub(pattern="^H",replacement="",x=colnames(MCMC.H))
    ## Remove trailing characters.
    hierarchical<-sub(pattern="\\[.*$",replacement="",x=hierarchical)
    ## Format values as numeric.
    hierarchical<-as.numeric(hierarchical)
    ## Get number of hierarchical variables.
    hierarchical<-max(hierarchical)
    ## Check number of hierarchical variables.
    if(J!=hierarchical){
      stop("Improper number of hierarchical variables.")
    }
    
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
    
    # If hierarchical variables are provided.
    if(J > 0){
      
      # Loop through each observation.
      for(i in 1:N){
        
        # Loop through each hierarchical variable.
        for(j in 1:J){
          
          # Retrieve hierarchical index.
          index<-H[i,j]
          
          # Create pattern for field names.
          fields<-paste0("^H",j,"\\[[[:digit:]]+,",index,"\\]$")
          
          # Retrieve hierarchical effect matrix.
          V<-MCMC.H[s,grepl(pattern=fields,x=colnames(MCMC.H)),
                    drop=FALSE]
          
          # Check for matching dimensions.
          if(ncol(V)!=D){
            stop("Dimension mismatch in hierarchical effects.")
          }
          
          # Apply hierarchical effects.
          eta[i,]<-eta[i,]+V
          
        }
        
      }
      
    }
    
    # Store linear combination.
    store[[s]]<-eta
    
  }
  
  # Create observation-wise storage list.
  P<-vector(mode="list",length=N)
  
  # Loop through each observation.
  for(i in 1:N){
    
    # Retrieve posterior predictions.
    P[[i]]<-t(sapply(X=store,FUN=function(x) x[i,]))
    
  }
  
  # Normalize posterior predictions with softmax.
  P<-sapply(X=P,FUN=softmax,simplify=FALSE)
  
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