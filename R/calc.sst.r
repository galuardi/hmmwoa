calc.sst <- function(tagdata, sst.dir, g, dateVec, raster = 'stack'){
  # compare tag sst data to oi sst map and calculate likelihoods
  
  #' @param: tagdata is variable containing tag-collected SST data
  #' @param: sst.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between tag-based sst and oi sst maps
  
  # get unique time points
  udates <- unique(sst$Date)
  
  for(i in 1:length(udates)){
    
    time <- udates[i]
    sst.i <- sst[which(sst$Date == time),]
    
    # open day's hycom data
    nc <- open.ncdf(paste(sst.dir, as.Date(time))) #add lat lon in filename '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'sst')
    lon <- get.var.ncdf(nc, 'lon')
    lat <- get.var.ncdf(nc, 'lat')
  
    # calculate sdx
    # ???
    
    # compare hycom to that day's tag-based ohc
    lik <- dnorm(sst, tag.sst, sdx) 
    
    if(i == 1){
      # result will be array of likelihood surfaces
      L.sst <- array(0, dim = c(dim(lik), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.sst[,,idx] = lik
    
    print(paste(time, ' finished.', sep=''))
    
  }
  
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.sst <- list(x = lon-360, y = lat, z = L.sst)
  ex <- extent(list.sst)
  L.ohc <- brick(list.sst$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.ohc <- flip(L.sst, direction = 'y')

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
