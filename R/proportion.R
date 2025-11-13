#' Grouped Proportion Plot
#'
#' @description Generates proportion plots for multiple groups.
#' @details Produces a pie-chart-like proportion plot with grouping structure. Each circle represents a group, and each sector represents a sample. When `a` = `TRUE` (the default), then the proportion of each sector filled with color represents the within-sample proportional abundance. Groups are sorted alphabetically and arranged from left to right and top to bottom. When `s` = `FALSE` (the default), then samples are sorted alphabetically and arranged in a clockwise orientation (from angle zero). When `s` = `TRUE`, then samples are sorted by decreasing proportional abundance. Samples are sorted independently for each group. This plot design is specialized for visualizing proportional abundance data.
#' @param x A list of vectors named `"g"`, `"s"`, and `"p"`. The elements of vector `"g"` (character or numeric) specify the group. The elements of vector `"s"` (character or numeric) specify the sample. The elements of vector `"p"` (numeric) specify the proportional abundance within sample `"s"`.
#' @param r Numeric scalar. Radius of plot circle (default = `1`).
#' @param b Numeric scalar. Plot radius buffer (proportion; default = `0.025`).
#' @param v Numeric scalar. Vertex count of plot circle (default = `1000`).
#' @param w Numeric scalar. Line width of outer circle (default = `1`).
#' @param f Numeric scalar. Line width of sectors as a proportion of `w` (default = `0.5`).
#' @param c Character string. Fill color of sector proportions (default = `"lightskyblue"`).
#' @param s Logical value. If `FALSE` (the default), sort samples alphabetically. If `TRUE`, sort samples by decreasing proportional abundance. Samples are sorted independently for each group.
#' @param a Logical value. If `FALSE`, proportional abundance is represented by the fraction of filled sector radius to outer sector radius. If `TRUE` (the default), proportional abundance is represented by the fraction of filled sector area to outer sector area.
#' @param m Numeric scalar. Maximum number of plot columns (default = `3`).
#' @param ... Additional arguments passed to [`title`][graphics::title()].
#' @returns No return value.
#' @seealso
#' [`singular.proportion`][singular.proportion()] for singular proportion plots. \cr \cr
#' [`detection`][detection()] for grouped detection plots.
#' @references A manuscript describing this plot design is in preparation.
#' @examplesIf interactive()
#' set.seed(1234)
#' n.groups<-6
#' n.samples<-6
#' data<-list(g=rep(x=LETTERS[1:n.groups],each=n.samples),
#'            s=rep(x=letters[1:n.samples],times=n.groups),
#'            p=stats::rbeta(n=n.groups*n.samples,
#'                           shape1=1,shape2=1))
#' proportion(x=data)
#' @export
proportion<-function(x,r=1,b=0.025,v=1e3,w=1,f=0.5,c="lightskyblue",s=FALSE,a=TRUE,m=3,...){
  
  # Check data input values.
  ## Ensure data object is a list.
  if(!is.list(x)) stop("Data must be a list.")
  ## Ensure data names include "g", "s", and "p".
  if(!all(c("g","s","p") %in% names(x))){
    stop("Data names must include: 'g', 's', 'p'.")
  }
  ## Subset data to vectors "g", "s", and "p".
  x<-x[c("g","s","p")]
  ## Ensure data contains vector elements.
  if(!all(sapply(X=x,FUN=is.vector))) stop("Non-vector elements in data.")
  ## Ensure group vector is numeric or character.
  if(!(is.numeric(x[["g"]]) | is.character(x[["g"]]))){
    stop("Group vector must be numeric or character.")
  }
  ## Ensure sample vector is numeric or character.
  if(!(is.numeric(x[["s"]]) | is.character(x[["s"]]))){
    stop("Sample vector must be numeric or character.")
  }
  ## Ensure proportions vector is numeric.
  if(!is.numeric(x[["p"]])) stop("Proportion vector must be numeric.")
  ## Ensure data vectors do not contain NAs.
  if(any(sapply(X=x,FUN=function(x) any(is.na(x))))) stop("Data contains NAs.")
  # Retrieve length of group vector.
  N<-length(x[["g"]])
  ## Ensure data vectors are the same length.
  if(!all(N==sapply(X=x,FUN=length))) stop("Length mismatch in data.")
  ## Ensure that proportions are between zero and one.
  if(!all((x[["p"]] >= 0) & (x[["p"]] <= 1))) stop("Invalid proportion values.")
  
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
  
  # Check input values for sort.
  ## Ensure sort vector is logical.
  if(!is.logical(s)) stop("Sort must be logical.")
  ## Ensure sort vector has length 1.
  if(length(s)!=1) stop("Sort must be of length 1.")
  ## Ensure sort value is not NA.
  if(is.na(s)) stop("Sort cannot be NA.")
  
  # Check input values for area.
  ## Ensure area vector is logical.
  if(!is.logical(a)) stop("Area must be logical.")
  ## Ensure area vector has length 1.
  if(length(a)!=1) stop("Area must be of length 1.")
  ## Ensure area value is not NA.
  if(is.na(a)) stop("Area cannot be NA.")
  
  # Check input values for column.
  ## Ensure column vector is numeric.
  if(!is.numeric(m)) stop("Column must be numeric.")
  ## Ensure column vector has length 1.
  if(length(m)!=1) stop("Column must be of length 1.")
  ## Ensure column value is not NA.
  if(is.na(m)) stop("Column cannot be NA.")
  
  # Sort group values.
  o<-sort(unique(x[["g"]]))
  
  # Retrieve number of groups.
  l<-length(o)
  
  # Update plot columns.
  if(l < m) m<-l
  
  # Derive plot rows.
  n<-ceiling(l/m)
  
  # Set plot array.
  graphics::par(mfrow=c(n,m))
  
  # Loop through each group.
  for(k in o){
    
    # Retrieve group indices.
    y<-as.numeric(x[["g"]]==k)
    
    # Subset data to group.
    u<-sapply(X=x,FUN=function(x) x[y==1],
              simplify=FALSE)
    
    # Generate singular plot for each group.
    singular.proportion(x=u,r=r,b=b,v=v,w=w,f=f,c=c,t=k,s=s,a=a,...)
    
  }
  
  # Reset plot array.
  graphics::par(mfrow=c(1,1))
  
}