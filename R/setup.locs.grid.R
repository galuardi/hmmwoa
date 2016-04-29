#' Setup the discrete spatial grid for the HMM
#' 
#' @param limits can take either a 'data.frame' of -Locations from WC psat tag
#'   or a pre-determined list of limits as \code{list(xmin, xmax, ymin, ymax)}
#' @param res character indicating resolution of grid. 'hycom' is .08 to match
#'   hycom reanalysis res. 'quarter' and 'one' are .25 and 1 deg, respectively.
#'   
#' @return a list
#' @export
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#'   
#' @examples
#' none
#' 
setup.locs.grid <- function(limits, res = 'quarter'){

  if(class(limits) == 'list'){
    # Find longitude extents
    il <- floor(limits[[1]])
    al <- ceiling(limits[[2]])
    lx <- 0.1 * (al - il)
    lonl <- il - lx
    lonu <- al + lx
    
    # Find latitude extents
    ila <- floor(limits[[3]])
    ala <- ceiling(limits[[4]])
    ly <- 0.1 * (ala - ila)
    latl <- ila - ly
    latu <- ala + ly
    
  } else if(class(limits) == 'data.frame'){

    # Find longitude extents
    il <- floor(min(locations$Longitude))
    al <- ceiling(max(locations$Longitude))
    lx <- 0.1 * (al - il)
    lonl <- il - lx
    lonu <- al + lx
    
    # Find latitude extents
    ila <- floor(min(locations$Latitude))
    ala <- ceiling(max(locations$Latitude))
    ly <- 0.1 * (ala - ila)
    latl <- ila - ly
    latu <- ala + ly
    
  } else{
    stop('Unable to determine format of input data.')
  }
  
  if(res == 'hycom'){
    lo.out <- 25/2 * (lonu - lonl)
    la.out <- 25/2 * (latu - latl) # 0.08 deg
  } else if(res == 'quarter'){
    lo.out <- 4 * (lonu - lonl)
    la.out <- 4 * (latu - latl) # .25 deg
  } else if(res == 'one'){
    lo.out <- 1 * (lonu - lonl)
    la.out <- 1 * (latu - latl) # 1 deg
  }
  
  # Create grid
  lo <- seq(lonl, lonu, length.out = lo.out)
  la <- seq(latl, latu, length.out = la.out)
  g <- meshgrid(lo, la)
  dlo <- lo[2] - lo[1]
  dla <- la[2] - la[1]
  
  list(lon = g$X, lat = g$Y, dlo = dlo, dla = dla)
}
