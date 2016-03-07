likint2 <- function(woa, woasd, minT, maxT){
  
  wlist = array(1e-6, dim=c(dim(woa)[1], dim(woa)[2], 2))
  wlist[,,1] = woa
  wlist[,,2] = woasd
  wlist[is.na(wlist)] = 1e-6
  as.matrix(aaply(wlist, 1:2, .fun = function(x) integrate(dnorm, lower = minT, upper = maxT , mean = x[1], sd = x[2])$value))
                  #.parallel=T, .paropts=list(.packages=c('stats','plyr'))))

}

#aaply(.parallel=T, .paropts=list(.packages=c('stats','plyr'), .verbose=T))