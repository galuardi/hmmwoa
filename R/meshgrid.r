
meshgrid <- function(x, y){
  # funtion from Pedersen et al 2011
  Y <- repmat(as.matrix(y), 1, length(x))
  X <- repmat(t(as.matrix(x)), length(y), 1)
  list(X = X, Y = Y)
}
