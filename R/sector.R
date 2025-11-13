#' Draw Sector Polygon
#'
#' @description Draws sector polygon.
#' @details Draws a sector polygon given a start angle, end angle, and circle radius. The sector is drawn about the origin (*i.e.*, x = 0, y = 0). Intended for use with [`template`][template()] to generate [`detection`][detection()] and [`proportion`][proportion] plots.
#' @param s Numeric scalar of start angle (degrees).
#' @param e Numeric scalar of end angle (degrees).
#' @param r Numeric scalar of circle radius.
#' @param v Numeric scalar of full-circle vertex count (default = `1000`).
#' @param ... Additional arguments passed to [`polgyon`][graphics::polygon()].
#' @returns No return value.
#' @seealso
#' [`circle`][circle()] for plotting circle polygons.
#' @examplesIf interactive()
#' template(l=1)
#' sector(s=0,e=45,r=1)
#' @export
sector<-function(s,e,r,v=1e3,...){
  
  # Check input values for start.
  ## Ensure start vector is numeric.
  if(!is.numeric(s)) stop("Start must be numeric.")
  ## Ensure start vector has length 1.
  if(length(s)!=1) stop("Start must be of length 1.")
  ## Ensure start value is not NA.
  if(is.na(s)) stop("Start cannot be NA.")
  
  # Check input values for end.
  ## Ensure end vector is numeric.
  if(!is.numeric(e)) stop("End must be numeric.")
  ## Ensure end vector has length 1.
  if(length(e)!=1) stop("End must be of length 1.")
  ## Ensure end value is not NA.
  if(is.na(e)) stop("End cannot be NA.")
  
  # Check input values for radius.
  ## Ensure radius vector is numeric.
  if(!is.numeric(r)) stop("Radius must be numeric.")
  ## Ensure radius vector has length 1.
  if(length(r)!=1) stop("Radius must be of length 1.")
  ## Ensure radius value is not NA.
  if(is.na(r)) stop("Radius cannot be NA.")
  
  # Check input values for vertex.
  ## Ensure vertex vector is numeric.
  if(!is.numeric(v)) stop("Vertex must be numeric.")
  ## Ensure vertex vector has length 1.
  if(length(v)!=1) stop("Vertex must be of length 1.")
  ## Ensure vertex value is not NA.
  if(is.na(v)) stop("Vertex cannot be NA.")
  
  # Generate angle vector.
  a<-seq(from=0,to=360,by=360/v)
  
  # Subset to applied range.
  a<-a[(a >= s) & (a <= e)]
  
  # Append start angle.
  if(a[1]!=s) a<-c(s,a)
  
  # Append end angle.
  if(a[length(a)]!=e) a<-c(a,e)
  
  # Generate coordinate matrix.
  m<-coordinates(a=a,r=r)
  
  # Append origin coordinates.
  m<-rbind(c(0,0),m,c(0,0))
  
  # Plot sector polygon.
  graphics::polygon(x=m,...)
  
}