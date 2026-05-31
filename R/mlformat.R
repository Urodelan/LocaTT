#' Formatting Multivariate Logistic Regression Models
#'
#' @description Format output of multivariate logistic regression models.
#' @details Formats output of a multivariate logistic regression model fit using the [`mlreg`][mlreg()] function. When `mode = "B"` (the default), returns regression coefficient estimates. When `mode = "R"`, returns residual correlation estimates (if [`mlreg`][mlreg()] was fit with `multivariate = TRUE`). Summarizes parameters by the quantiles of their posterior distributions, with a point estimate at the 50th percentile (*i.e.*, the posterior median). Lower and upper limits are defined by the credible interval argument. At the default `ci = 0.95`, returns 95% credible intervals. When a credible interval does not overlap zero, the point estimate is appended with an asterisk.
#' @param fit A `stanfit` object returned from the [`mlreg`][mlreg()] function. The fitted multivariate logistic regression model.
#' @param mode Character scalar. Specifies which parameter set to summarize. When `mode = "B"` (the default), generates summaries of regression coefficient estimates. When `mode = "R"`, generates summaries of residual correlation estimates (if [`mlreg`][mlreg()] was fit with `multivariate = TRUE`).
#' @param ci Numeric scalar. Defines the credible interval for parameter summaries. When `ci = 0.95` (the default), generates 95% credible intervals for parameter estimates. `ci` must be between `0` and `1`.
#' @param digits Numeric scalar. Positive integer value specifying the number of decimal places to which results will be rounded. The default is `3`.
#' @param names.x Character vector (optional). If provided, then supplies the names of predictor variables in the returned matrix. Names should match those of the predictor variables used to fit the [`mlreg`][mlreg()] model.
#' @param names.y Character vector (optional). If provided, then supplies the names of response variables in the returned matrix. Names should match those of the response variables used to fit the [`mlreg`][mlreg()] model.
#' @returns Numeric matrix of posterior summaries.
#' @seealso
#' [`mlreg`][mlreg()] for fitting multivariate logistic regression models. \cr \cr
#' [`mlcoef`][mlcoef()] for extracting regression coefficients from multivariate logistic regression models. \cr \cr
#' [`mlcor`][mlcor()] for extracting residual correlations from multivariate logistic regression models.
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
#' # Retrieve fitted regression model.
#' fit<-data$fit
#' 
#' # Retrieve predictor matrix.
#' X<-data$X
#' 
#' # Retrieve response matrix.
#' Y<-data$Y
#' 
#' # Extract regression coefficients.
#' B<-mlformat(fit=fit,mode="B",
#'             names.x=colnames(X),
#'             names.y=colnames(Y))
#' 
#' # Display regression coefficients.
#' print(B,quote=FALSE,right=TRUE)
#' 
#' # Extract residual correlations.
#' R<-mlformat(fit=fit,mode="R",
#'             names.y=colnames(Y))
#' 
#' # Display residual correlations.
#' print(R,quote=FALSE,right=TRUE)
#' @export
mlformat<-function(fit,mode="B",ci=0.95,digits=3,names.x,names.y){
  
  # Check that model is a stanfit object.
  if(!inherits(x=fit,what="stanfit")) stop("Model must be a 'stanfit' object.")
  
  # Check that model is a mlreg fit.
  if(fit@model_name!="mlreg") stop("Model must be a 'mlreg' fit.")
  
  # Check mode input data.
  ## Check that mode is a vector.
  if(!is.vector(mode)) stop("mode must be a vector.")
  ## Check that mode is not a list.
  if(is.list(mode)) stop("mode cannot be a list.")
  ## Check that mode has length 1.
  if(length(mode)!=1) stop("mode must have length 1.")
  ## Check that mode is not NA.
  if(is.na(mode)) stop("mode cannot be NA.")
  ## Check that mode is character.
  if(!is.character(mode)) stop("mode must be character.")
  ## Check that ci is B or R.
  if(!(mode %in% c("B","R"))) stop("mode must be 'B' or 'R'.")
  
  # Check ci input data.
  ## Check that ci is a vector.
  if(!is.vector(ci)) stop("ci must be a vector.")
  ## Check that ci is not a list.
  if(is.list(ci)) stop("ci cannot be a list.")
  ## Check that ci has length 1.
  if(length(ci)!=1) stop("ci must have length 1.")
  ## Check that ci is not NA.
  if(is.na(ci)) stop("ci cannot be NA.")
  ## Check that ci is numeric.
  if(!is.numeric(ci)) stop("ci must be numeric.")
  ## Check that ci is between 0 and 1.
  if((ci < 0) | (ci > 1)) stop("ci must be between 0 and 1.")
  
  # Check digits input data.
  ## Check that digits is a vector.
  if(!is.vector(digits)) stop("digits must be a vector.")
  ## Check that digits is not a list.
  if(is.list(digits)) stop("digits cannot be a list.")
  ## Check that digits has length 1.
  if(length(digits)!=1) stop("digits must have length 1.")
  ## Check that digits is not NA.
  if(is.na(digits)) stop("digits cannot be NA.")
  ## Check that digits is numeric.
  if(!is.numeric(digits)) stop("digits must be numeric integer.")
  ## Check that digits is integer.
  if(digits %% 1 != 0) stop("digits must be numeric integer.")
  ## Check that digits is positive.
  if(digits < 0) stop("digits must be positive.")
  
  # Check names.x input data.
  if(!missing(names.x)){
    ## Check that names.x is a vector.
    if(!is.vector(names.x)) stop("names.x must be a vector.")
    ## Check that names.x is not a list.
    if(is.list(names.x)) stop("names.x cannot be a list.")
    ## Check that names.x is character.
    if(!is.character(names.x)) stop("names.x must be character.")
    ## Issue warning if mode is R.
    if(mode=="R") warning("names.x ignored when mode = 'R'.")
  }
  
  # Check names.y input data.
  if(!missing(names.y)){
    ## Check that names.y is a vector.
    if(!is.vector(names.y)) stop("names.y must be a vector.")
    ## Check that names.y is not a list.
    if(is.list(names.y)) stop("names.y cannot be a list.")
    ## Check that names.y is character.
    if(!is.character(names.y)) stop("names.y must be character.")
  }
  
  # Define retrieval function.
  retrieval<-function(x,digits){
    ## Retrieve lower quantile.
    q.l<-x[1]
    ## Retrieve middle quantile.
    q.m<-x[2]
    ## Retrieve upper quantile.
    q.u<-x[3]
    ## Check for significance.
    sig<-as.numeric((q.l > 0) | (q.u < 0))
    ## Format lower quantile.
    q.l<-formatC(x=q.l,digits=digits,format="f")
    ## Format middle quantile.
    q.m<-formatC(x=q.m,digits=digits,format="f")
    ## Format upper quantile.
    q.u<-formatC(x=q.u,digits=digits,format="f")
    ## Append significance identifier.
    if(sig) q.m<-paste0(q.m,"*")
    ## Collect information.
    out<-paste0(q.m," (",q.l," to ",q.u,")")
    ## Return information.
    return(out)
  }
  
  # Compute bounds.
  ## Lower bound.
  l<-(1-ci)/2
  ## Upper bound.
  u<-1-l
  
  # Retrieve summaries.
  if(mode=="B"){ ## When mode is B.
    ## Retrieve regression coefficients.
    B<-mlcoef(fit=fit,probs=c(l,0.5,u))
    ## Extract output by array dimension names.
    out<-apply(X=B,MARGIN=c("X","Y"),FUN=retrieval,digits=digits)
  }else{ ## When mode is R.
    ## Retrieve residual correlations.
    R<-mlcor(fit=fit,probs=c(l,0.5,u))
    ## Get indices of dimensions labeled as Y.
    indices<-which(names(dimnames(R))=="Y")
    ## Extract output by array dimension indices.
    out<-apply(X=R,MARGIN=indices,FUN=retrieval,digits=digits)
  }
  
  # Finish formatting.
  ## Clear dimension names.
  names(dimnames(out))<-NULL
  ## Transpose summaries.
  out<-t(out)
  ## Apply provided names.
  if(mode=="B"){ ## When mode is B.
    ## If names.y is provided.
    if(!missing(names.y)){
      ## Verify number of records.
      if(nrow(out)!=length(names.y)) stop("names.y has improper length.")
      ## Apply row names.
      row.names(out)<-names.y
    }
    ## If names.x is provided.
    if(!missing(names.x)){
      ## Verify number of fields.
      if(ncol(out)!=length(names.x)) stop("names.x has improper length.")
      ## Apply field names.
      colnames(out)<-names.x
    }
  }else{ ## When mode is R.
    ## If names.y is provided.
    if(!missing(names.y)){
      ## Verify number of records.
      if(nrow(out)!=length(names.y)) stop("names.y has improper length.")
      ## Verify number of fields.
      if(ncol(out)!=length(names.y)) stop("names.y has improper length.")
      ## Apply row names.
      row.names(out)<-names.y
      ## Apply field names.
      colnames(out)<-names.y
    }
  }
  
  # Return summaries.
  return(out)
  
}