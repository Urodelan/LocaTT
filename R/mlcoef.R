#' Coefficients of Multivariate Logistic Regression Models
#'
#' @description Extract regression coefficients from multivariate logistic regression models.
#' @details Extracts regression coefficient estimates from a multivariate logistic regression model fit using the [`mlreg`][mlreg()] function. Summarizes estimates by the quantiles of their posterior distributions, and returns summaries in a 3-dimensional array. The dimensions of the 3D array represent the posterior quantiles, the predictor variables, and the response variables, respectively. Values at `probs = 0.025` and `0.975` comprise 95% credible intervals. Values at `probs = 0.25` and `0.75` comprise 50% credible intervals, and values at `probs = 0.5` represent point estimates.
#' @param fit A `stanfit` object returned from the [`mlreg`][mlreg()] function. The fitted multivariate logistic regression model.
#' @param probs Numeric vector of probabilities. Passed to the `probs` argument of the [`stats::quantile`][stats::quantile()] function. The length of `probs` defines the length of the first dimension of the returned 3-dimensional array. By default, `probs = c(0.025, 0.25, 0.5, 0.75, 0.975)`. See details for the interpretation of values at each probability.
#' @param dimnames List (optional). If provided, then names within the returned 3-dimensional array will receive these values. Passed to the `dimnames` argument of the [`array`][array()] function. If omitted, then generic names will be provided to the returned 3D array. See the `dimnames` argument of the [`array`][array()] function for details.
#' @returns Numeric 3-dimensional array of regression coefficient posterior quantiles.
#' @seealso
#' [`mlreg`][mlreg()] for fitting multivariate logistic regression models. \cr \cr
#' [`mlcor`][mlcor()] for extracting residual correlations from multivariate logistic regression models. \cr \cr
#' [`mlformat`][mlformat()] for formatting output of multivariate logistic regression models.
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
#' # Extract regression coefficients.
#' out<-mlcoef(fit=data$fit)
#' @export
mlcoef<-function(fit,probs=c(0.025,0.25,0.5,0.75,0.975),dimnames){
  
  # Check that model is a stanfit object.
  if(!inherits(x=fit,what="stanfit")) stop("Model must be a 'stanfit' object.")
  
  # Check that model is a mlreg fit.
  if(fit@model_name!="mlreg") stop("Model must be a 'mlreg' fit.")
  
  # Check probs input data.
  ## Check that probs is a vector.
  if(!is.vector(probs)) stop("probs must be a vector.")
  ## Check that probs are not a list.
  if(is.list(probs)) stop("probs cannot be a list.")
  ## Check that probs does not contain NAs.
  if(any(is.na(probs))) stop("probs cannot contain NAs.")
  ## Check that probs is numeric.
  if(!is.numeric(probs)) stop("probs must be numeric.")
  ## Check that probs is between 0 and 1.
  if(any((probs < 0) | (probs > 1))) stop("probs must be between 0 and 1.")
  
  # Retrieve posterior samples.
  MCMC<-as.matrix(fit)
  
  # Subset to regression coefficients.
  MCMC.B<-MCMC[,grepl(pattern="^B\\[",x=colnames(MCMC)),drop=FALSE]
  
  # Retrieve number of predictor variables.
  ## Removing leading characters.
  K<-sub(pattern="^B\\[[[:digit:]]+,",replacement="",x=colnames(MCMC.B))
  ## Remove trailing characters.
  K<-sub(pattern="\\]$",replacement="",x=K)
  ## Format values as numeric.
  K<-as.numeric(K)
  ## Get number of predictor variables.
  K<-max(K)
  
  # Retrieve number of response variables.
  ## Removing leading characters.
  D<-sub(pattern="^B\\[",replacement="",x=colnames(MCMC.B))
  ## Remove trailing characters.
  D<-sub(pattern="\\,[[:digit:]]+]$",replacement="",x=D)
  ## Format values as numeric.
  D<-as.numeric(D)
  ## Get number of response variables.
  D<-max(D)
  
  # Get length of probabilities.
  N<-length(probs)
  
  # Compute quantiles.
  Q<-t(apply(X=MCMC.B,MARGIN=2,FUN=stats::quantile,probs=probs))
  
  # Re-transpose for single probability.
  if(N==1) Q<-t(Q)
  
  # If dimension names are missing.
  if(missing(dimnames)){
    
    # Initialize digits storage vector.
    digits<-c()
    
    # Loop through each probability.
    for(i in 1:N){
      
      # Retrieve probability.
      prob<-probs[i]
      
      # If probability is an integer.
      if(prob %% 1 == 0){
        
        # Then the probability has no digits.
        digit<-0
        
      }else{ # If probability is not an integer.
        
        # Retrieve digits for the probability.
        digit<-nchar(sub(pattern="^.*\\.",replacement="",x=prob))
        
      }
      
      # Store probability digits.
      digits<-c(digits,digit)
      
    }
    
    # Get maximum number of probability digits.
    digits<-max(digits)
    
    # Generate dimension names.
    ## 1st dimension: Probability.
    names.N<-paste0("Q",formatC(x=probs,digits=digits,format="f"))
    ## 2nd dimension: Predictor.
    names.K<-paste0("X",formatC(x=1:K,width=nchar(K),format="d",flag="0"))
    ## 3rd dimension: Response.
    names.D<-paste0("Y",formatC(x=1:D,width=nchar(D),format="d",flag="0"))
    
    # Collect dimension names.
    dimnames<-list(names.N,names.K,names.D)
    
    # Provide names for dimnames.
    names(dimnames)<-c("Q","X","Y")
    
  }
  
  # Initialize three-dimensional array.
  A<-array(dim=c(N,K,D),dimnames=dimnames)
  
  # Loop through each probability.
  for(i in 1:N){
    
    # Retrieve regression coefficient matrix.
    A[i,,]<-matrix(data=Q[,i],nrow=K,ncol=D,byrow=TRUE)
    
  }
  
  # Return regression coefficient estimates.
  return(A)
  
}