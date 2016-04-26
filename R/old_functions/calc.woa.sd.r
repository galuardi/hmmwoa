# calculate sd on a woa grid

#' Title calc.woa.sd
#' Calcualtes standard deviation according to surrounding cells
#' @param woa 
#' @references LeBris, A., Fréchet, A., Galbraith, P.S., and Wroblewski, J.S. 2013. Evidence for alternative migratory behaviours in the northern Gulf of St Lawrence population of Atlantic cod (Gadus morhua L.). ICES J. Mar. Sci. doi:10.1093/icesjms/fst068.

#' @return a grid of sd values
#' @export
#'
#' @examples
calc.woa.sd <- function(woa, extent = 3){
  r = flip(raster(t(dat)),2)
  f1 = focal(r, w = matrix(1/(extent*extent),nrow = extent, ncol = extent), fun = sd)
  f1
}


