#' Initiate Template Plot
#'
#' @description Initiates a blank template plot.
#' @details Initiates a blank template plot with limits `l` and buffer `b` about the origin (*i.e.*, x = 0, y = 0). `l` is used for axis limits in both the negative and positive directions. `b` extends the limits beyond `l` by a fixed proportion (*i.e.*, `l` * (1 + `b`)). Intended for use with [`circle`][circle()] and [`sector`][sector()].
#' @param l Numeric scalar of axis limits (applies to both axes).
#' @param b Numeric scalar to extend axis limits (see Details; default = `0.025`).
#' @returns No return value.
#' @seealso
#' [`circle`][circle()] for plotting circle polygons. \cr \cr
#' [`sector`][sector()] for plotting sector polygons.
#' @examplesIf interactive()
#' template(l=1)
#' circle(r=1)
#' @export
template<-function(l,b=0.025){
  
  # Check input values for limit.
  ## Ensure limit vector is numeric.
  if(!is.numeric(l)) stop("Limit must be numeric.")
  ## Ensure limit vector has length 1.
  if(length(l)!=1) stop("Limit must be of length 1.")
  ## Ensure limit value is not NA.
  if(is.na(l)) stop("Limit cannot be NA.")
  
  # Check input values for buffer.
  ## Ensure buffer vector is numeric.
  if(!is.numeric(b)) stop("Buffer must be numeric.")
  ## Ensure buffer vector has length 1.
  if(length(b)!=1) stop("Buffer must be of length 1.")
  ## Ensure buffer value is not NA.
  if(is.na(b)) stop("Buffer cannot be NA.")
  
  # Generate empty template plot.
  graphics::plot.default(x=0,y=0,type="p",col=NA, # Transparent point at origin.
                         xlim=c(-l,l)*(1+b),ylim=c(-l,l)*(1+b), # Buffer extent.
                         asp=1,xaxs="i",yaxs="i",bty="n", # Format axes.
                         axes=FALSE,ann=FALSE) # Remove annotations.
  
}