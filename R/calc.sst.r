#' Calculate SST-based likelihood
#' 
#' \code{calc.sst} compares tag SST to remotely sensed SST and calculates
#' likelihoods
#' 
#' @param tag.sst variable containing tag-collected SST data
#' @param sst.dir local directory where remote sensing SST downloads are stored
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#'   
#' @return likelihood is raster brick of likelihood surfaces representing matches
#'   between tag-based sst and remotely sensed sst maps
#' @export
#' @seealso \code{\link{calc.ohc}}
#' @examples
#' # see example script blue256_example.r

calc.sst <- function(tag.sst, ptt, sst.dir, dateVec){
  
  start.t <- Sys.time()
  
  dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
  tag.sst[,12] <- as.Date(dts)
  names(tag.sst)[12] <- 'dts'
  by_dte <- dplyr::group_by(tag.sst, tag.sst$dts)  # group by unique DAILY time points
  tag.sst <- data.frame(dplyr::summarise(by_dte, min(by_dte$Temperature), max(by_dte$Temperature)))
  colnames(tag.sst) <- list('date', 'minT', 'maxT')
  
  T <- length(tag.sst[,1])
  
  print(paste('Starting iterations through time ', '...'))
  
  for(i in 1:T){
    
    time <- tag.sst$date[i]
    sst.i <- c(tag.sst$minT[i] * .95, tag.sst$maxT[i] * 1.05) # sensor error
    
    # open day's sst data
    nc <- RNetCDF::open.nc(paste(sst.dir, ptt, '_', as.Date(time), '.nc', sep='')) #add lat lon in filename '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, 'sst') # for OI SST
    
    # calc sd of SST
    # focal calc on mean temp and write to sd var
    r = raster::flip(raster::raster(t(dat)), 2)
    sdx = raster::focal(r, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) stats::sd(x, na.rm = T))
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))

    # compare sst to that day's tag-based ohc
    lik.sst <- likint3(dat, sdx, sst.i[1], sst.i[2])

    if(i == 1){
      lon <- RNetCDF::var.get.nc(nc, 'longitude')
      lat <- RNetCDF::var.get.nc(nc, 'latitude')
      # result will be array of likelihood surfaces
      L.sst <- array(0, dim = c(dim(lik.sst), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.sst[,,idx] = (lik.sst / max(lik.sst, na.rm=T)) - 0.2
    
  }
    
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.sst <- list(x = lon, y = lat, z = L.sst)
  ex <- raster::extent(list.sst)
  L.sst <- raster::brick(list.sst$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.sst <- raster::flip(L.sst, direction = 'y')
    
  L.sst[L.sst < 0] <- 0
  
  # return sst likelihood surfaces
  return(L.sst)
    
}
  