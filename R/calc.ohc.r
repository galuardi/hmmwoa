calc.ohc <- function(tagdata, isotherm = '', ohc.dir, dateVec, raster = T){
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
  
  ohcVec <- rep(0, length.out = length(udates))
  
  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  for(i in 1:length(udates)){
    
      time <- udates[i]
      pdt.i <- pdt[which(pdt$Date == time),]
      
      # open day's hycom data
      nc <- open.ncdf(paste(ohc.dir, 'Lyd_', as.Date(time), '.nc', sep=''))
      depth <- get.var.ncdf(nc, 'depth')
      lon <- get.var.ncdf(nc, 'lon')
      lat <- get.var.ncdf(nc, 'lat')
      
      # isotherm is minimum temperature recorded for that time point
      if(iso.def == FALSE) isotherm <- min(pdt.i$MinTemp, na.rm = T)
      
      # perform tag data integration
      tag <- approx(pdt.i$Depth, pdt.i$MidTemp, xout = depth)
      tag <- tag$y - isotherm
      tag.ohc <- cp * rho * sum(tag, na.rm = T) / 10000
      
      # store tag ohc
      ohcVec[i] <- tag.ohc
      
  }
  
  for(i in 1:length(udates)){

    # define time based on tag data
    time <- udates[i]
    
    # open day's hycom data
    nc <- open.ncdf(paste(ohc.dir, 'Lyd_', as.Date(time), '.nc', sep=''))
    dat <- get.var.ncdf(nc, 'water_temp')
    #depth <- get.var.ncdf(nc, 'depth')
    #lon <- get.var.ncdf(nc, 'lon')
    #lat <- get.var.ncdf(nc, 'lat')
    
    pdt.i <- pdt[which(pdt$Date == time),]
    
    # calculate daily isotherm based on tag data
    if(iso.def == FALSE) isotherm <- min(pdt.i$MinTemp, na.rm = T)
    
    dat[dat<isotherm] <- NA
    
    # Perform hycom integration
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat, 1:2, sum, na.rm = T) / 10000 
    
    if(i == 1){
      sdx <- sd(ohcVec[c(1,2)])
    } else{
      sdx <- sd(ohcVec[c((i - 1), i, (i + 1))])
    }
    
    # compare hycom to that day's tag-based ohc
    #lik.dt <- matrix(dtnorm(ohc, tag.ohc, sdx, 0, 150), dim(ohc)[1], dim(ohc)[2])
    lik <- dnorm(ohc, ohcVec[i], sdx) 
    #lik <- (lik / max(lik, na.rm = T)) - .05 # normalize
    #print(paste(max(lik), time))
    
    if(i == 1){
      # result should be array of likelihood surfaces
      L.ohc <- array(0, dim = c(dim(lik), length(dateVec)))
    }
    
    idx <- which(dateVec == as.Date(time))
    L.ohc[,,idx] = lik
    
    print(paste(time, ' finished.', sep=''))
    
  }
  
  if(raster){
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    list.ohc <- list(x = lon-360, y = lat, z = L.ohc)
    ex <- extent(list.ohc)
    L.ohc <- brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
    L.ohc <- flip(L.ohc, direction = 'y')
    L.ohc <- stack(L.ohc)
  }
  
  print(class(L.ohc))
  # return ohc likelihood surfaces
  return(L.ohc)
  
}
