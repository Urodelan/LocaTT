#' Generate Circular Coordinates
#'
#' @description Generates coordinates along a circular path.
#' @details Calculates xy coordinates along a circular path given a vector of angles and a specified radius. The center of the circle is at the origin (*i.e.*, x = 0, y = 0). This function is helpful for calculating the vertices of circle and sector polygons.
#' @param a Numeric vector of angles (degrees).
#' @param r Numeric scalar of the circle radius.
#' @returns A numeric matrix of xy coordinates.
#' @seealso
#' [`circle`][circle()] for plotting circle polygons. \cr \cr
#' [`sector`][sector()] for plotting sector polygons.
#' @examples
#' coordinates(a=c(90,180,270,360),r=1)
#' @export
coordinates<-function(a,r){
  
  # Check input values for angle.
  ## Ensure angle vector is numeric.
  if(!is.numeric(a)) stop("Angle must be numeric.")
  ## Ensure angle vector does not contain NAs.
  if(any(is.na(a))) stop("Angle cannot contain NAs.")
  
  # Check input values for radius.
  ## Ensure radius vector is numeric.
  if(!is.numeric(r)) stop("Radius must be numeric.")
  ## Ensure radius vector has length 1.
  if(length(r)!=1) stop("Radius must be of length 1.")
  ## Ensure radius value is not NA.
  if(is.na(r)) stop("Radius cannot be NA.")
  
  # Define pi.
  pi<-acos(-1)
  
  # Define internal functions.
  ## Define x-coordinate function.
  x<-function(a,r) r*sin(pi/180*a)
  ## Define y-coordinate function.
  y<-function(a,r) r*cos(pi/180*a)
  ## Define xy-coordinate function.
  xy<-function(a,r) c(x=x(a,r),y=y(a,r))
  
  # Generate coordinates.
  m<-t(sapply(X=a,FUN=xy,r=r))
  
  # Return coordinates.
  return(m)
  
}