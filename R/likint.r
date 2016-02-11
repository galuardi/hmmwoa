likint <- function(woa, minT, maxT, sdT){
  matrix(vapply(woa, FUN = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x, sd = sdT)$value, FUN.VALUE = matrix(0,1,1)), dim(woa)[1], dim(woa)[2])
}
