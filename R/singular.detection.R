#' Singular Detection Plot
#'
#' @description Generates a detection plot for a singular group.
#' @details Produces a pie-chart-like detection plot without grouping structure. Each sector represents a sample, and each sub-sector represents a replicate. Filled replicates represent detections. Samples are sorted alphabetically and arranged in a clockwise orientation (from angle zero). This plot design is specialized for visualizing binary detection data.
#' @param x A list of vectors named `"s"`, `"r"`, and `"d"`. The elements of vector `"s"` (character or numeric) specify the sample. The elements of vector `"r"` (numeric) specify the number of replicates within sample `"s"`. The elements of vector `"r"` (numeric) specify the number of replicates within sample `"s"` with detections.
#' @param r Numeric scalar. Radius of plot circle (default = `1`).
#' @param b Numeric scalar. Plot radius buffer (proportion; default = `0.025`).
#' @param v Numeric scalar. Vertex count of plot circle (default = `1000`).
#' @param w Numeric scalar. Line width of outer circle (default = `1`).
#' @param f Numeric scalar. Line width of sectors as a proportion of `w` (default = `0.5`).
#' @param c Character string. Fill color of sub-sector detections (default = `"lightskyblue"`).
#' @param t Character string. Plot title (default = `""`).
#' @param ... Additional arguments passed to [`title`][graphics::title()].
#' @returns No return value.
#' @seealso
#' [`detection`][detection()] for grouped detection plots. \cr \cr
#' [`proportion`][proportion()] for grouped proportion plots.
#' @references A manuscript describing this plot design is in preparation.
#' @examplesIf interactive()
#' set.seed(1234)
#' n.samples<-6
#' n.replicates<-3
#' data<-list(s=letters[1:n.samples],
#'            r=rep(x=n.replicates,times=n.samples),
#'            d=sample(x=0:n.replicates,size=n.samples,
#'                     replace=TRUE))
#' singular.detection(x=data)
#' @export
singular.detection<-function(x,r=1,b=0.025,v=1e3,w=1,f=0.5,c="lightskyblue",t="",...){
  
  # Check data input values.
  ## Ensure data object is a list.
  if(!is.list(x)) stop("Data must be a list.")
  ## Ensure data names include "s", "r", and "d".
  if(!all(c("s","r","d") %in% names(x))){
    stop("Data names must include: 's', 'r', 'd'.")
  }
  ## Subset data to vectors "s", "r", and "d".
  x<-x[c("s","r","d")]
  ## Ensure data contains vector elements.
  if(!all(sapply(X=x,FUN=is.vector))) stop("Non-vector elements in data.")
  ## Ensure sample vector is numeric or character.
  if(!(is.numeric(x[["s"]]) | is.character(x[["s"]]))){
    stop("Sample vector must be numeric or character.")
  }
  ## Ensure data vectors are numeric.
  if(!all(sapply(X=x[c("r","d")],FUN=is.numeric))) stop("Non-numeric elements in data.")
  ## Ensure data vectors do not contain NAs.
  if(any(sapply(X=x,FUN=function(x) any(is.na(x))))) stop("Data contains NAs.")
  # Retrieve number of samples.
  N<-length(x[["s"]])
  ## Ensure data vectors are the same length.
  if(!all(N==sapply(X=x,FUN=length))) stop("Length mismatch in data.")
  ## Ensure that detections do not exceed replicates.
  if(!all(x[["d"]] <= x[["r"]])) stop("Detections exceed replicates.")
  ## Retrieve sample order.
  o<-order(x[["s"]])
  ## Sort data by sample identifier.
  x<-sapply(X=x,FUN=function(x) x[o],
            simplify=FALSE)
  
  # Check input values for radius.
  ## Ensure radius vector is numeric.
  if(!is.numeric(r)) stop("Radius must be numeric.")
  ## Ensure radius vector has length 1.
  if(length(r)!=1) stop("Radius must be of length 1.")
  ## Ensure radius value is not NA.
  if(is.na(r)) stop("Radius cannot be NA.")
  
  # Check input values for buffer.
  ## Ensure buffer vector is numeric.
  if(!is.numeric(b)) stop("Buffer must be numeric.")
  ## Ensure buffer vector has length 1.
  if(length(b)!=1) stop("Buffer must be of length 1.")
  ## Ensure buffer value is not NA.
  if(is.na(b)) stop("Buffer cannot be NA.")
  
  # Check input values for width.
  ## Ensure width vector is numeric.
  if(!is.numeric(w)) stop("Width must be numeric.")
  ## Ensure width vector has length 1.
  if(length(w)!=1) stop("Width must be of length 1.")
  ## Ensure width value is not NA.
  if(is.na(w)) stop("Width cannot be NA.")
  
  # Check input values for fraction.
  ## Ensure fraction vector is numeric.
  if(!is.numeric(f)) stop("Fraction must be numeric.")
  ## Ensure fraction vector has length 1.
  if(length(f)!=1) stop("Fraction must be of length 1.")
  ## Ensure fraction value is not NA.
  if(is.na(f)) stop("Fraction cannot be NA.")
  ## Ensure fraction value is between 0 and 1.
  if((f <= 0) | (f > 1)) stop("Fraction must be between 0 and 1.")
  
  # Check input values for vertex.
  ## Ensure vertex vector is numeric.
  if(!is.numeric(v)) stop("Vertex must be numeric.")
  ## Ensure vertex vector has length 1.
  if(length(v)!=1) stop("Vertex must be of length 1.")
  ## Ensure vertex value is not NA.
  if(is.na(v)) stop("Vertex cannot be NA.")
  
  # Check input values for color.
  ## Ensure color vector is character.
  if(!is.character(c)) stop("Color must be character.")
  ## Ensure color vector has length 1.
  if(length(c)!=1) stop("Color must be of length 1.")
  ## Ensure color value is not NA.
  if(is.na(c)) stop("Color cannot be NA.")
  
  # Check input values for title.
  ## Ensure title vector is character.
  if(!is.character(t)) stop("Title must be character.")
  ## Ensure title vector has length 1.
  if(length(t)!=1) stop("Title must be of length 1.")
  ## Ensure title value is not NA.
  if(is.na(t)) stop("Title cannot be NA.")
  
  # Initialize plot.
  template(l=r,b=b)
  
  # Generate angle increments.
  u<-seq(from=0,to=360,by=360/N)
  
  # Loop through each sample.
  for(i in 1:N){
    
    # Retrieve number of replicates.
    s.r<-x[["r"]][i]
    
    # Retrieve number of detections.
    s.d<-x[["d"]][i]
    
    # Loop through each replicate.
    for(j in s.r:1){
      
      # Set color given detection.
      color<-ifelse(j <= s.d,c,"white")
      
      # Plot sector polygon.
      sector(s=u[i],e=u[i+1],r=r/s.r*j,
             v=v,col=color,lwd=w*f)
      
    }
    
  }
  
  # Plot circle polygon.
  circle(r=r,v=v,lwd=w)
  
  # Produce plot title.
  graphics::title(main=t,font.main=1,...)
  
}