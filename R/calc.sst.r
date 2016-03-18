calc.sst <- function(tagdata, sst.dir, g = g, dateVec, raster = 'stack'){
  # compare tag sst data to oi sst map and calculate likelihoods
  
  #' @param: tagdata is variable containing tag-collected SST data
  #' @param: sst.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between tag-based sst and oi sst maps
  require('dplyr')
  
  dts <- as.POSIXct(tagdata$Date, format = findDateFormat(tagdata$Date))
  
  tagdata[,12] <- as.Date(dts)
  by_dte <- group_by(tagdata, V12)  # group by unique DAILY time points
  tagdata <- data.frame(summarise(by_dte, min(Temperature), max(Temperature)))
  colnames(tagdata) <- list('date','minT','maxT')
  
  T <- length(tagdata[,1])
  
  for(i in 1:T){
    
    time <- tagdata$date[i]
    sst.i <- c(tagdata$minT[i] * .99, tagdata$maxT[i] * 1.01) # sensor error
    
    # open day's sst data
    nc <- open.ncdf(paste(sst.dir, ptt, '_', as.Date(time), '.nc', sep='')) #add lat lon in filename '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'sst') # for OI SST
    
    # calc sd of SST
    # focal calc on mean temp and write to sd var
    t = Sys.time()
    r = flip(raster(t(dat)),2)
    sdx = focal(r, w=matrix(1,nrow=3,ncol=3), fun=function(x) sd(x, na.rm = T))
    sdx = t(as.matrix(flip(sdx,2)))
    print(paste('finishing sd for ', time,'. Section took ', Sys.time() - t))
    
    # compare sst to that day's tag-based ohc
    t = Sys.time()
    lik.sst <- likint3(dat, sdx, sst.i[1], sst.i[2])
    print(paste('finishing lik.sst for ', time,'. Section took ', Sys.time() - t))
    
    if(i == 1){
      lon <- get.var.ncdf(nc, 'longitude')
      lat <- get.var.ncdf(nc, 'latitude')
      # result will be array of likelihood surfaces
      L.sst <- array(0, dim = c(dim(lik.sst), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.sst[,,idx] = lik.sst
  }
    
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    list.sst <- list(x = lon, y = lat, z = L.sst)
    ex <- extent(list.sst)
    L.sst <- brick(list.sst$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    L.sst <- flip(L.sst, direction = 'y')
    
    # make L.sst match resolution/extent of g
    row <- dim(g$lon)[1]
    col <- dim(g$lon)[2]
    ex <- extent(c(min(g$lon[1,]), max(g$lon[1,]), min(g$lat[,1]), max(g$lat[,1])))
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    rasMatch <- raster(ex, nrows=row, ncols=col, crs = crs)
    L.sst <- spatial_sync_raster(L.sst, rasMatch)
    
    if(raster == 'brick'){
      s <- L.sst
    } else if(raster == 'stack'){
      s <- stack(L.sst)
    } else if(raster == 'array'){
      s <- raster::as.array(L.sst, transpose = T)
    }
    
    print(class(L.sst))
    # return sst likelihood surfaces
    return(L.sst)
    
}
  