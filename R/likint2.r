likint2 <- function(woa, woasd, minT, maxT, intLib = 'pracma'){
  
  wlist = array(1e-6, dim=c(dim(woa)[1], dim(woa)[2], 2))
  wlist[,,1] = woa
  wlist[,,2] = woasd
  wlist[is.na(wlist)] = 1e-6
  if(intLib=='pracma'){
    as.matrix(aaply(wlist, 1:2, .fun = function(x) pracma::integral(dnorm, xmin = minT, xmax = maxT , mean = x[1], sd = x[2])))
  } else if(intLib=='stats'){
    as.matrix(aaply(wlist, 1:2, .fun = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x[1], sd = x[2])$value))
  } else{
    stop('No integration library specified.')
  }
                
}