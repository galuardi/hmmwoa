
repmat <- function(X, m, n){
  # function from Pedersen et al 2011
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}