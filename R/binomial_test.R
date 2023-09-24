#' Binomial Test
#'
#' @description Performs binomial tests.
#' @details Calls on the [`pbinom`][stats::pbinom()] function in the [stats][stats::stats] package to perform vectorized binomial tests. Arguments are recycled as in [`pbinom`][stats::pbinom()]. Only one-sided tests are supported, and only p-values are returned.
#' @param k A numeric vector of the number of successes.
#' @param n A numeric vector of the number of trials.
#' @param p A numeric vector of the hypothesized probabilities of success.
#' @param alternative A string specifying the alternative hypothesis. Must be `"less"` or `"greater"` (the default).
#' @returns A numeric vector of p-values from the binomial tests.
#' @examples
#' binomial_test(k=c(5,1,7,4),
#'               n=c(10,3,15,5),
#'               p=c(0.2,0.1,0.5,0.6),
#'               alternative="greater")
#' @export
binomial_test<-function(k,n,p,alternative="greater"){
  
  # Throw an error if the number of successes are not numeric.
  if(!is.numeric(k)) stop("The number of successes (k) must be numeric.")

  # Throw an error if the number of trials are not numeric.
  if(!is.numeric(n)) stop("The number of trials (n) must be numeric.")

  # Throw an error if the hypothesized probabilities of success are not numeric.
  if(!is.numeric(p)) stop("The hypothesized probabilities of success (p) must be numeric.")
  
  # Throw an error if alternative is not a character string.
  if(!is.character(alternative)) stop("Alternative must be a character string.")
  # Throw an error if the length of the alternative character string is not one.
  if(length(alternative)!=1) stop("Alternative must be a character string of length one.")
  # Throw an error if alternative is not "less" or "greater".
  if(!(alternative %in% c("less","greater"))) stop("Alternative must be 'less' or 'greater'.")
  
  # If the alternative hypothesis is "less".
  if(alternative=="less"){
    
    # Calculate Pr(X <= k) if p is true.
    p.values<-stats::pbinom(q=k,size=n,prob=p,lower.tail=TRUE,log.p=FALSE)
    
  } else { # If the alternative hypothesis is "greater".
    
    # Calculate Pr(X >= k) if p is true. Note that one is subtracted from k
    # because stats::pbinom with lower.tail=FALSE calculates Pr(X > q).
    p.values<-stats::pbinom(q=k-1,size=n,prob=p,lower.tail=FALSE,log.p=FALSE)
    
  }
  
  # Return the p-values.
  return(p.values)
  
}