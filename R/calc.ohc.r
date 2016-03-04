calc.ohc <- function(tagdata, isotherm = '', ohc.dir, g, dateVec, raster = 'stack'){
  # compare tag data to ohc map and calculate likelihoods
  
  #' @param: tagdata is variable containing tag-collected PDT data
  #' @param: isotherm is default '' in which isotherm is calculated
  #' on the fly based on daily tag data. Otherwise, numeric isotherm
  #' constraint can be specified (e.g. 20).
  #' @param: ohc.dir is local directory where get.hycom downloads are
  #' stored.
  #' @return: likelihood is array of likelihood surfaces representing
  #' matches between tag-based ohc and hycom ohc maps
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  # get unique time points
  udates <- unique(pdt$Date)
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  for(i in 1:length(udates)){
    
    time <- udates[i]
    pdt.i <- pdt[which(pdt$Date == time),]
    
    # open day's hycom data
    nc <- open.ncdf(paste(ohc.dir, 'Lyd_', as.Date(time), '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'water_temp')
    depth <- get.var.ncdf(nc, 'depth')
    lon <- get.var.ncdf(nc, 'lon')
    lat <- get.var.ncdf(nc, 'lat')
    
    if(i==1){
      f.arr <- array(NA, dim=c(length(lon),length(lat),length(udates)))
    }
    
    # isotherm is minimum temperature recorded for that time point
    if(iso.def == FALSE) isotherm <- min(pdt.i$MinTemp, na.rm = T)
    
    # perform tag data integration
    tag <- approx(pdt.i$Depth, pdt.i$MidTemp, xout = depth)
    tag <- tag$y - isotherm
    tag.ohc <- cp * rho * sum(tag, na.rm = T) / 10000
    
    # calc ohc for min/max temps for each day to calc sdx for dnorm
    minTag <- approx(pdt.i$Depth, pdt.i$MinTemp, xout = depth)
    minTag <- minTag$y - isotherm
    minT.ohc <- cp * rho * sum(minTag, na.rm = T) / 10000
    
    maxTag <- approx(pdt.i$Depth, pdt.i$MaxTemp, xout = depth)
    maxTag <- maxTag$y - isotherm
    maxT.ohc <- cp * rho * sum(maxTag, na.rm = T) / 10000
    
    #sdx <- sd(c(minT.ohc, maxT.ohc))
    #print(sdx)
    
    # Perform hycom integration
    dat[dat<isotherm] <- NA
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat, 1:2, sum, na.rm = T) / 10000 
    
    # calc sd of OHC
    # focal calc on mean temp and write to sd var
    list.r <- list(x = lon, y = lat, z = ohc)
    ex <- extent(list.r)
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    r <- raster(t(list.r$z), xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
    r <- flip(r, direction = 'y')
    # focal matrix and calculation
    w = matrix(1/9, nrow = 3, ncol = 3)
    f <- focal(r, w, function(x) sd(x))
    
    # put results in an array
    f.arr[,,i] <- t(as.matrix(flip(f,direction='y')))
    sdx <- f.arr[,,i]
    # compare hycom to that day's tag-based ohc
    lik <- dnorm(tag.ohc, mean=ohc, sd=sdx) 
    
    if(i == 1){
      # result will be array of likelihood surfaces
      L.ohc <- array(0, dim = c(dim(lik), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.ohc[,,idx] = lik
    
    print(paste(time, ' finished.', sep=''))
    
  }
  
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.ohc <- list(x = lon-360, y = lat, z = L.ohc)
  ex <- extent(list.ohc)
  L.ohc <- brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.ohc <- flip(L.ohc, direction = 'y')
  s <- stack(L.ohc)
  return(s)
}
  # make L.pdt match resolution/extent of g
  #row <- dim(g$lon)[1]
  #col <- dim(g$lon)[2]
  #ex <- extent(c(min(g$lon[1,]), max(g$lon[1,]), min(g$lat[,1]), max(g$lat[,1])))
  #crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  #rasMatch <- raster(ex, nrows=row, ncols=col, crs = crs)
  #L.ohc <- spatial_sync_raster(L.ohc, rasMatch)
  
  #if(raster == 'brick'){
   # s <- L.ohc
  #} else if(raster == 'stack'){
  #  s <- stack(L.ohc)
  #} else if(raster == 'array'){
  #  s <- raster::as.array(L.ohc, transpose = T)
  #}
  
#  print(class(L.ohc))
  # return ohc likelihood surfaces
 # return(L.ohc)
  
#}
