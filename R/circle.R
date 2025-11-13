#' Draw Circle Polygon
#'
#' @description Draws circle polygon.
#' @details Draws a circle polygon of a given radius. The circle is drawn about the origin (*i.e.*, x = 0, y = 0). Intended for use with [`template`][template()] to generate [`detection`][detection()] and [`proportion`][proportion] plots.
#' @param r Numeric scalar of circle radius.
#' @param v Numeric scalar of vertex count (default = `1000`).
#' @param ... Additional arguments passed to [`polgyon`][graphics::polygon()].
#' @returns No return value.
#' @seealso
#' [`sector`][sector()] for plotting sector polygons.
#' @examplesIf interactive()
#' template(l=1)
#' circle(r=1)
#' @export
circle<-function(r,v=1e3,...){
  
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
  
  # Generate coordinate matrix.
  m<-coordinates(a=a,r=r)
  
  # Plot circle polygon.
  graphics::polygon(x=m,...)
  
}