#' Repeat your matrix?
#' 
#' @param x ...
#' @param m ...
#' @param n ...
#' @return matrix ...
#'   
#' @examples
#' x <- c(1, 2, 3)
#' repmat(x)
#' 
#' @references
#' Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#' Estimating animal behavior and residency from movement data. Oikos
#' 120:1281–1290. doi: 10.1111/j.1600-0706.2011.19044.x

repmat <- function(X, m, n){
  # function from Pedersen et al 2011
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}