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
#' @seealso \code{\link{calc.pdt.int}}, \code{\link{calc.ohc}}
#' @examples
#' # see example script blue256_example.r

calc.sst <- function(tag.sst, sst.dir, dateVec){
  
  start.t <- Sys.time()
  
  dts <- as.POSIXct(tag.sst$Date, format = findDateFormat(tag.sst$Date))
  tag.sst[,12] <- as.Date(dts)
  by_dte <- dplyr::group_by(tag.sst, V12)  # group by unique DAILY time points
  tag.sst <- data.frame(dplyr::summarise(by_dte, min(Temperature), max(Temperature)))
  colnames(tag.sst) <- list('date', 'minT', 'maxT')
  
  T <- length(tag.sst[,1])
  
  for(i in 1:T){
    
    time <- tag.sst$date[i]
    sst.i <- c(tag.sst$minT[i] * .99, tag.sst$maxT[i] * 1.01) # sensor error
    
    # open day's sst data
    nc <- ncdf::open.ncdf(paste(sst.dir, ptt, '_', as.Date(time), '.nc', sep='')) #add lat lon in filename '.nc', sep=''))
    dat <- ncdf::get.var.ncdf(nc, 'sst') # for OI SST
    
    # calc sd of SST
    # focal calc on mean temp and write to sd var
    t = Sys.time()
    r = raster::flip(raster::raster(t(dat)), 2)
    sdx = raster::focal(r, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) sd(x, na.rm = T))
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
    print(paste('finishing sd for ', time,'. Section took ', Sys.time() - t))
    
    # compare sst to that day's tag-based ohc
    t = Sys.time()
    lik.sst <- likint3(dat, sdx, sst.i[1], sst.i[2])
    print(paste('finishing lik.sst for ', time,'. Section took ', Sys.time() - t))
    
    if(i == 1){
      lon <- ncdf::get.var.ncdf(nc, 'longitude')
      lat <- ncdf::get.var.ncdf(nc, 'latitude')
      # result will be array of likelihood surfaces
      L.sst <- array(0, dim = c(dim(lik.sst), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.sst[,,idx] = lik.sst
  }
    
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    list.sst <- list(x = lon, y = lat, z = L.sst)
    ex <- raster::extent(list.sst)
    L.sst <- raster::brick(list.sst$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    L.sst <- raster::flip(L.sst, direction = 'y')
    
    print(paste('SST calculations took ', Sys.time() - start.t, '.'))
    
    # return sst likelihood surfaces
    return(L.sst)
    
}
  