#' Calculate density function and integrate between limits
#' 
#' #' \code{likint3} calculates density function for a normal distribution and integrates between limits
#'  
#' @param w is an array of grid values like sea surface temperature
#' @param wsd is an array containing the sd of the grid values, usually from \code{raster::focal}
#' @param minT is an integer representing the lower limit of the tag-measured variable (e.g. SST)
#' @param maxT is an integer representing the upper limit of the tag-measured variable
#' 
#' @return an array of dim(w) that represents the likelihood of the tag-measured variable as compared to the input grid
#'   

likint3 <- function(w, wsd, minT, maxT){
  midT = (maxT + minT) / 2
  Tsd = (maxT - minT) / 4
  widx = w >= minT & w <= maxT & !is.na(w)
  wdf = data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
  wdf$wsd[is.na(wdf$wsd)] = 0
  # wint = apply(wdf, 1, function(x) pracma::integral(dnorm, minT, maxT, mean = x[1], sd = x[2]))
  wint = apply(wdf, 1, function(x) stats::integrate(stats::dnorm, x[1]-x[2], x[1]+x[2], mean = midT, sd = Tsd * 2)$value) 
  w = w * 0
  w[widx] = wint
  w
} 